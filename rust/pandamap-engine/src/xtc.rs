//! Pure-Rust GROMACS XTC trajectory reader and writer.
//!
//! Implements the exact XDR + `xdrfile` 3dfcoord compressed coordinate codec
//! (magic 1995 standard XTC and magic 2023 XTC64) with zero C/GROMACS
//! dependencies. Coordinates are returned in Angstrom (`nm * 10.0`).

use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;

const XTC_MAGIC: i32 = 1995;
const XTC64_MAGIC: i32 = 2023;
const FIRSTIDX: usize = 9;
const LASTIDX: usize = 73;

#[rustfmt::skip]
const MAGICINTS: [u32; 73] = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 10, 12, 16, 20, 25, 32, 40, 50, 64,
    80, 101, 128, 161, 203, 256, 322, 406, 512, 645, 812, 1024, 1290,
    1625, 2048, 2580, 3250, 4096, 5060, 6501, 8192, 10321, 13003,
    16384, 20642, 26007, 32768, 41285, 52015, 65536, 82570, 104031,
    131072, 165140, 208063, 262144, 330280, 416127, 524287, 660561,
    832255, 1048576, 1321122, 1664510, 2097152, 2642245, 3329021,
    4194304, 5284491, 6658042, 8388607, 10568983, 13316085, 16777216,
];

/// Metadata + raw compressed payload of one XTC frame (allows parallel decoding).
#[derive(Debug, Clone)]
pub struct RawXtcFrame {
    pub frame_idx: usize,
    pub natoms: usize,
    pub step: i32,
    pub time_ps: f32,
    pub box_angstrom: [[f32; 3]; 3],
    pub payload: XtcPayload,
}

#[derive(Debug, Clone)]
pub enum XtcPayload {
    Uncompressed(Vec<f32>), // raw nm floats, length 3 * natoms
    Compressed {
        precision: f32,
        minint: [i32; 3],
        maxint: [i32; 3],
        smallidx: usize,
        buf: Vec<u8>,
    },
}

/// Decoded frame coordinates in Angstrom (`[x, y, z]` per atom).
#[derive(Debug, Clone)]
pub struct XtcFrame {
    pub frame_idx: usize,
    pub step: i32,
    pub time_ps: f32,
    pub box_angstrom: [[f32; 3]; 3],
    pub coords_angstrom: Vec<[f32; 3]>,
}

impl RawXtcFrame {
    /// Decompress this frame into Angstrom coordinates (`nm * 10.0`).
    pub fn decode(&self) -> io::Result<XtcFrame> {
        let coords_angstrom = match &self.payload {
            XtcPayload::Uncompressed(raw) => {
                let mut out = Vec::with_capacity(self.natoms);
                for c in raw.chunks_exact(3) {
                    out.push([c[0] * 10.0, c[1] * 10.0, c[2] * 10.0]);
                }
                out
            }
            XtcPayload::Compressed {
                precision,
                minint,
                maxint,
                smallidx,
                buf,
            } => decompress_3dfcoord(self.natoms, *precision, *minint, *maxint, *smallidx, buf)?,
        };
        Ok(XtcFrame {
            frame_idx: self.frame_idx,
            step: self.step,
            time_ps: self.time_ps,
            box_angstrom: self.box_angstrom,
            coords_angstrom,
        })
    }
}

pub struct XtcReader<R: Read + Seek> {
    reader: R,
    next_frame_idx: usize,
}

impl XtcReader<BufReader<File>> {
    pub fn open<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::open(path)?;
        Ok(Self::new(BufReader::with_capacity(256 * 1024, file)))
    }
}

impl<R: Read + Seek> XtcReader<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            next_frame_idx: 0,
        }
    }

    /// Read the next raw frame from the XTC stream. Returns `Ok(None)` at EOF.
    pub fn read_next_raw(&mut self) -> io::Result<Option<RawXtcFrame>> {
        let magic = match read_i32_be_opt(&mut self.reader)? {
            Some(m) => m,
            None => return Ok(None),
        };
        if magic != XTC_MAGIC && magic != XTC64_MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Invalid XTC magic number: {} (expected 1995 or 2023)", magic),
            ));
        }

        let natoms_i32 = read_i32_be(&mut self.reader)?;
        if natoms_i32 < 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Negative XTC atom count: {}", natoms_i32),
            ));
        }
        let natoms = natoms_i32 as usize;
        let step = read_i32_be(&mut self.reader)?;
        let time_ps = read_f32_be(&mut self.reader)?;

        let mut box_angstrom = [[0.0f32; 3]; 3];
        for row in &mut box_angstrom {
            for val in row.iter_mut() {
                *val = read_f32_be(&mut self.reader)? * 10.0;
            }
        }

        let lsize = read_i32_be(&mut self.reader)? as usize;
        if lsize != natoms {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("XTC lsize {} != natoms {}", lsize, natoms),
            ));
        }

        let payload = if natoms <= 9 {
            let mut raw = vec![0.0f32; natoms * 3];
            for v in &mut raw {
                *v = read_f32_be(&mut self.reader)?;
            }
            XtcPayload::Uncompressed(raw)
        } else {
            let precision = read_f32_be(&mut self.reader)?;
            let minint = [
                read_i32_be(&mut self.reader)?,
                read_i32_be(&mut self.reader)?,
                read_i32_be(&mut self.reader)?,
            ];
            let maxint = [
                read_i32_be(&mut self.reader)?,
                read_i32_be(&mut self.reader)?,
                read_i32_be(&mut self.reader)?,
            ];
            let smallidx = read_i32_be(&mut self.reader)? as usize;
            let nbytes: usize = if magic == XTC64_MAGIC {
                read_i64_be(&mut self.reader)? as usize
            } else {
                read_i32_be(&mut self.reader)? as usize
            };
            let mut buf = vec![0u8; nbytes];
            self.reader.read_exact(&mut buf)?;
            let pad = (4 - (nbytes % 4)) % 4;
            if pad > 0 {
                let mut skip = [0u8; 4];
                self.reader.read_exact(&mut skip[..pad])?;
            }
            XtcPayload::Compressed {
                precision,
                minint,
                maxint,
                smallidx,
                buf,
            }
        };

        let frame_idx = self.next_frame_idx;
        self.next_frame_idx += 1;
        Ok(Some(RawXtcFrame {
            frame_idx,
            natoms,
            step,
            time_ps,
            box_angstrom,
            payload,
        }))
    }

    /// Skip the next frame without allocating or decompressing coordinates.
    #[allow(dead_code)]
    pub fn skip_next(&mut self) -> io::Result<Option<(usize, f32)>> {
        let magic = match read_i32_be_opt(&mut self.reader)? {
            Some(m) => m,
            None => return Ok(None),
        };
        if magic != XTC_MAGIC && magic != XTC64_MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Invalid XTC magic number: {}", magic),
            ));
        }
        let natoms = read_i32_be(&mut self.reader)? as usize;
        let _step = read_i32_be(&mut self.reader)?;
        let time_ps = read_f32_be(&mut self.reader)?;
        self.reader.seek(SeekFrom::Current(9 * 4 + 4))?; // box + lsize
        if natoms <= 9 {
            self.reader.seek(SeekFrom::Current((natoms * 3 * 4) as i64))?;
        } else {
            self.reader.seek(SeekFrom::Current(4 + 3 * 4 + 3 * 4 + 4))?;
            let nbytes: i64 = if magic == XTC64_MAGIC {
                read_i64_be(&mut self.reader)?
            } else {
                read_i32_be(&mut self.reader)? as i64
            };
            let padded = (nbytes + 3) & !3;
            self.reader.seek(SeekFrom::Current(padded))?;
        }
        let idx = self.next_frame_idx;
        self.next_frame_idx += 1;
        Ok(Some((idx, time_ps)))
    }
}

#[inline]
fn read_i32_be_opt<R: Read>(r: &mut R) -> io::Result<Option<i32>> {
    let mut b = [0u8; 4];
    let mut total = 0;
    while total < 4 {
        let n = r.read(&mut b[total..])?;
        if n == 0 {
            if total == 0 {
                return Ok(None);
            }
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "Truncated XTC header",
            ));
        }
        total += n;
    }
    Ok(Some(i32::from_be_bytes(b)))
}

#[inline]
fn read_i32_be<R: Read>(r: &mut R) -> io::Result<i32> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)?;
    Ok(i32::from_be_bytes(b))
}

#[inline]
fn read_i64_be<R: Read>(r: &mut R) -> io::Result<i64> {
    let mut b = [0u8; 8];
    r.read_exact(&mut b)?;
    Ok(i64::from_be_bytes(b))
}

#[inline]
fn read_f32_be<R: Read>(r: &mut R) -> io::Result<f32> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)?;
    Ok(f32::from_be_bytes(b))
}

#[inline]
fn sizeofint(size: u32) -> usize {
    if size == 0 {
        return 0;
    }
    (32 - size.leading_zeros()) as usize
}

#[inline]
fn sizeofints3(sizes: [u32; 3]) -> usize {
    let prod: u128 = (sizes[0] as u128) * (sizes[1] as u128) * (sizes[2] as u128);
    if prod == 0 {
        return 0;
    }
    (128 - prod.leading_zeros()) as usize
}

struct BitReader<'a> {
    buf: &'a [u8],
    cnt: usize,
    lastbits: u32,
    lastbyte: u64,
}

impl<'a> BitReader<'a> {
    #[inline]
    fn new(buf: &'a [u8]) -> Self {
        Self {
            buf,
            cnt: 0,
            lastbits: 0,
            lastbyte: 0,
        }
    }

    #[inline]
    fn receivebits(&mut self, mut nbits: usize) -> u32 {
        let mut num: u32 = 0;
        while nbits >= 8 {
            let b = *self.buf.get(self.cnt).unwrap_or(&0) as u64;
            self.cnt += 1;
            self.lastbyte = (self.lastbyte << 8) | b;
            num = (num << 8) | (((self.lastbyte >> self.lastbits) & 0xff) as u32);
            nbits -= 8;
        }
        if nbits > 0 {
            let nb = nbits as u32;
            if self.lastbits < nb {
                let b = *self.buf.get(self.cnt).unwrap_or(&0) as u64;
                self.cnt += 1;
                self.lastbyte = (self.lastbyte << 8) | b;
                self.lastbits += 8;
            }
            self.lastbits -= nb;
            let mask = (1u64 << nb) - 1;
            num = (num << nb) | (((self.lastbyte >> self.lastbits) & mask) as u32);
        }
        num
    }

    #[inline]
    fn receiveints3(&mut self, mut nbits: usize, sizes: [u32; 3]) -> [i32; 3] {
        let mut val: u128 = 0;
        let mut shift: u32 = 0;
        while nbits > 8 {
            let b = self.receivebits(8) as u128;
            val |= b << shift;
            shift += 8;
            nbits -= 8;
        }
        if nbits > 0 {
            let b = self.receivebits(nbits) as u128;
            val |= b << shift;
        }
        let s2 = sizes[2] as u128;
        let s1 = sizes[1] as u128;
        let n2 = (val % s2) as i32;
        val /= s2;
        let n1 = (val % s1) as i32;
        val /= s1;
        let n0 = val as i32;
        [n0, n1, n2]
    }
}

fn decompress_3dfcoord(
    natoms: usize,
    precision: f32,
    minint: [i32; 3],
    maxint: [i32; 3],
    mut smallidx: usize,
    buf: &[u8],
) -> io::Result<Vec<[f32; 3]>> {
    if precision == 0.0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "XTC frame has zero precision",
        ));
    }
    if smallidx >= LASTIDX {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Invalid XTC smallidx {}", smallidx),
        ));
    }

    let sizeint = [
        (maxint[0] as i64 - minint[0] as i64 + 1) as u32,
        (maxint[1] as i64 - minint[1] as i64 + 1) as u32,
        (maxint[2] as i64 - minint[2] as i64 + 1) as u32,
    ];
    let (bitsize, bitsizeint) = if (sizeint[0] | sizeint[1] | sizeint[2]) > 0x00ff_ffff {
        (
            0usize,
            [
                sizeofint(sizeint[0]),
                sizeofint(sizeint[1]),
                sizeofint(sizeint[2]),
            ],
        )
    } else {
        (sizeofints3(sizeint), [0, 0, 0])
    };

    let mut smallnum = (MAGICINTS[smallidx] / 2) as i32;
    let mut smaller = if smallidx > FIRSTIDX {
        (MAGICINTS[smallidx - 1] / 2) as i32
    } else {
        0
    };
    let mut sizesmall = [MAGICINTS[smallidx]; 3];

    // Convert directly to Angstrom: (int / precision) * 10.0
    let scale = 10.0f32 / precision;
    let mut out = Vec::with_capacity(natoms);
    let mut br = BitReader::new(buf);

    let mut i = 0usize;
    let mut run: usize = 0; // Persists across loop iterations when flag == 0 (exact xdrfile spec!)
    while i < natoms {
        let mut thiscoord = if bitsize == 0 {
            [
                br.receivebits(bitsizeint[0]) as i32,
                br.receivebits(bitsizeint[1]) as i32,
                br.receivebits(bitsizeint[2]) as i32,
            ]
        } else {
            br.receiveints3(bitsize, sizeint)
        };
        thiscoord[0] += minint[0];
        thiscoord[1] += minint[1];
        thiscoord[2] += minint[2];

        i += 1;
        let flag = br.receivebits(1);
        let mut is_smaller: i32 = 0;
        if flag == 1 {
            let raw_run = br.receivebits(5) as usize;
            let rem = raw_run % 3;
            run = raw_run - rem;
            is_smaller = rem as i32 - 1;
        }

        if run > 0 {
            let mut prevcoord = thiscoord;
            let mut k = 0usize;
            while k < run {
                let mut delta = br.receiveints3(smallidx, sizesmall);
                i += 1;
                delta[0] += prevcoord[0] - smallnum;
                delta[1] += prevcoord[1] - smallnum;
                delta[2] += prevcoord[2] - smallnum;
                if k == 0 {
                    // Interchange first and second atom (exact xdrfile water optimization)
                    let tmp = delta;
                    delta = prevcoord;
                    prevcoord = tmp;
                    out.push([
                        prevcoord[0] as f32 * scale,
                        prevcoord[1] as f32 * scale,
                        prevcoord[2] as f32 * scale,
                    ]);
                } else {
                    prevcoord = delta;
                }
                out.push([
                    delta[0] as f32 * scale,
                    delta[1] as f32 * scale,
                    delta[2] as f32 * scale,
                ]);
                k += 3;
            }
        } else {
            out.push([
                thiscoord[0] as f32 * scale,
                thiscoord[1] as f32 * scale,
                thiscoord[2] as f32 * scale,
            ]);
        }

        smallidx = (smallidx as i32 + is_smaller).clamp(FIRSTIDX as i32, (LASTIDX - 1) as i32) as usize;
        if is_smaller < 0 {
            smallnum = smaller;
            smaller = if smallidx > FIRSTIDX {
                (MAGICINTS[smallidx - 1] / 2) as i32
            } else {
                0
            };
        } else if is_smaller > 0 {
            smaller = smallnum;
            smallnum = (MAGICINTS[smallidx] / 2) as i32;
        }
        sizesmall = [MAGICINTS[smallidx]; 3];
    }

    if out.len() > natoms {
        out.truncate(natoms);
    }
    Ok(out)
}

// =========================================================================
// Spec-compliant GROMACS XTC Writer (useful for tests, benchmarks, & export)
// =========================================================================

struct BitWriter {
    buf: Vec<u8>,
    lastbits: u32,
    lastbyte: u64,
}

impl BitWriter {
    fn new() -> Self {
        Self {
            buf: Vec::new(),
            lastbits: 0,
            lastbyte: 0,
        }
    }

    fn sendbits(&mut self, mut nbits: usize, num: u32) {
        let mut shift = nbits;
        while nbits >= 8 {
            shift -= 8;
            let byte = ((num >> shift) & 0xff) as u64;
            self.lastbyte = (self.lastbyte << 8) | byte;
            self.buf.push(((self.lastbyte >> self.lastbits) & 0xff) as u8);
            nbits -= 8;
        }
        if nbits > 0 {
            let nb = nbits as u32;
            let mask = (1u64 << nb) - 1;
            let bits = (num as u64) & mask;
            self.lastbyte = (self.lastbyte << nb) | bits;
            self.lastbits += nb;
            if self.lastbits >= 8 {
                self.lastbits -= 8;
                self.buf.push(((self.lastbyte >> self.lastbits) & 0xff) as u8);
            }
        }
    }

    fn sendints3(&mut self, mut nbits: usize, sizes: [u32; 3], nums: [u32; 3]) {
        let mut val: u128 = nums[0] as u128;
        val = val * (sizes[1] as u128) + (nums[1] as u128);
        val = val * (sizes[2] as u128) + (nums[2] as u128);
        while nbits > 8 {
            self.sendbits(8, (val & 0xff) as u32);
            val >>= 8;
            nbits -= 8;
        }
        if nbits > 0 {
            self.sendbits(nbits, (val & ((1u128 << nbits) - 1)) as u32);
        }
    }

    fn finish(mut self) -> Vec<u8> {
        if self.lastbits > 0 {
            let pad = 8 - self.lastbits;
            self.buf.push(((self.lastbyte << pad) & 0xff) as u8);
        }
        self.buf
    }
}

pub struct XtcWriter<W: Write> {
    writer: W,
}

impl XtcWriter<BufWriter<File>> {
    pub fn create<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::create(path)?;
        Ok(Self {
            writer: BufWriter::with_capacity(256 * 1024, file),
        })
    }
}

impl<W: Write> XtcWriter<W> {
    pub fn new(writer: W) -> Self {
        Self { writer }
    }

    /// Write one frame; `coords_angstrom` and `box_angstrom` are converted to nm.
    pub fn write_frame(
        &mut self,
        step: i32,
        time_ps: f32,
        box_angstrom: [[f32; 3]; 3],
        coords_angstrom: &[[f32; 3]],
        precision: f32,
    ) -> io::Result<()> {
        let natoms = coords_angstrom.len();
        self.writer.write_all(&XTC_MAGIC.to_be_bytes())?;
        self.writer.write_all(&(natoms as i32).to_be_bytes())?;
        self.writer.write_all(&step.to_be_bytes())?;
        self.writer.write_all(&time_ps.to_be_bytes())?;
        for row in &box_angstrom {
            for &v in row {
                self.writer.write_all(&(v * 0.1).to_be_bytes())?;
            }
        }
        self.writer.write_all(&(natoms as i32).to_be_bytes())?;

        if natoms <= 9 {
            for c in coords_angstrom {
                self.writer.write_all(&(c[0] * 0.1).to_be_bytes())?;
                self.writer.write_all(&(c[1] * 0.1).to_be_bytes())?;
                self.writer.write_all(&(c[2] * 0.1).to_be_bytes())?;
            }
            return Ok(());
        }

        let scale = precision * 0.1; // Angstrom -> quantized nm integer
        let mut qcoords = Vec::with_capacity(natoms);
        let mut minint = [i32::MAX; 3];
        let mut maxint = [i32::MIN; 3];
        for c in coords_angstrom {
            let q = [
                (c[0] * scale).round() as i32,
                (c[1] * scale).round() as i32,
                (c[2] * scale).round() as i32,
            ];
            for d in 0..3 {
                if q[d] < minint[d] {
                    minint[d] = q[d];
                }
                if q[d] > maxint[d] {
                    maxint[d] = q[d];
                }
            }
            qcoords.push(q);
        }

        let sizeint = [
            (maxint[0] as i64 - minint[0] as i64 + 1) as u32,
            (maxint[1] as i64 - minint[1] as i64 + 1) as u32,
            (maxint[2] as i64 - minint[2] as i64 + 1) as u32,
        ];
        let (bitsize, bitsizeint) = if (sizeint[0] | sizeint[1] | sizeint[2]) > 0x00ff_ffff {
            (
                0usize,
                [
                    sizeofint(sizeint[0]),
                    sizeofint(sizeint[1]),
                    sizeofint(sizeint[2]),
                ],
            )
        } else {
            (sizeofints3(sizeint), [0, 0, 0])
        };

        // Write standalone atoms: first atom sets run=0 via flag=1 (raw_run=1 -> rem=1 -> run=0, is_smaller=0);
        // subsequent atoms use flag=0 (reusing run=0, saving 5 bits/atom).
        let smallidx: usize = FIRSTIDX;
        let mut bw = BitWriter::new();
        for (idx, q) in qcoords.iter().enumerate() {
            let rel = [
                (q[0] - minint[0]) as u32,
                (q[1] - minint[1]) as u32,
                (q[2] - minint[2]) as u32,
            ];
            if bitsize == 0 {
                bw.sendbits(bitsizeint[0], rel[0]);
                bw.sendbits(bitsizeint[1], rel[1]);
                bw.sendbits(bitsizeint[2], rel[2]);
            } else {
                bw.sendints3(bitsize, sizeint, rel);
            }
            if idx == 0 {
                bw.sendbits(1, 0); // initial run is 0, so flag=0 keeps run=0, is_smaller=0
            } else {
                bw.sendbits(1, 0);
            }
        }
        let buf = bw.finish();

        self.writer.write_all(&precision.to_be_bytes())?;
        for &m in &minint {
            self.writer.write_all(&m.to_be_bytes())?;
        }
        for &m in &maxint {
            self.writer.write_all(&m.to_be_bytes())?;
        }
        self.writer.write_all(&(smallidx as i32).to_be_bytes())?;
        self.writer.write_all(&(buf.len() as i32).to_be_bytes())?;
        self.writer.write_all(&buf)?;
        let pad = (4 - (buf.len() % 4)) % 4;
        if pad > 0 {
            self.writer.write_all(&[0u8; 4][..pad])?;
        }
        Ok(())
    }

    pub fn flush(&mut self) -> io::Result<()> {
        self.writer.flush()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn xtc_roundtrip_compressed_and_uncompressed() {
        // Test > 9 atoms (compressed 3dfcoord)
        let mut coords = Vec::new();
        for i in 0..25 {
            coords.push([10.0 + i as f32 * 0.75, -5.25 + i as f32 * 0.3, 2.5 + i as f32 * 0.1]);
        }
        let box_a = [[50.0, 0.0, 0.0], [0.0, 50.0, 0.0], [0.0, 0.0, 50.0]];
        let mut bytes = Vec::new();
        {
            let mut w = XtcWriter::new(&mut bytes);
            w.write_frame(0, 0.0, box_a, &coords, 1000.0).unwrap();
            let shifted: Vec<[f32; 3]> = coords.iter().map(|c| [c[0] + 1.2, c[1], c[2]]).collect();
            w.write_frame(5000, 1000.0, box_a, &shifted, 1000.0).unwrap();
            w.flush().unwrap();
        }
        let mut r = XtcReader::new(Cursor::new(bytes));
        let f0 = r.read_next_raw().unwrap().unwrap().decode().unwrap();
        let f1 = r.read_next_raw().unwrap().unwrap().decode().unwrap();
        assert!(r.read_next_raw().unwrap().is_none());
        assert_eq!(f0.step, 0);
        assert!((f1.time_ps - 1000.0).abs() < 1e-4);
        for (a, b) in coords.iter().zip(f0.coords_angstrom.iter()) {
            assert!((a[0] - b[0]).abs() < 0.02);
            assert!((a[1] - b[1]).abs() < 0.02);
            assert!((a[2] - b[2]).abs() < 0.02);
        }
        assert!((f1.coords_angstrom[0][0] - (coords[0][0] + 1.2)).abs() < 0.02);
    }
}
