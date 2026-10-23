"""Regression for SIGPROC extraction when optional nsamples is present/absent."""
import json
import struct
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

EXTRACTOR = Path(__file__).resolve().parents[1] / "scripts" / "extract_sigproc_windows.py"


def field(f, s):
    b=s.encode("utf-8")
    f.write(struct.pack("<i",len(b)))
    f.write(b)


def make_fil(path, with_nsamples):
    samples=41
    with path.open("wb") as f:
        field(f,"HEADER_START")
        for key,value,fmt in [
            ("source_name","TEST",None),
            ("tstart",60000.12345,"<d"),
            ("tsamp",0.001,"<d"),
            ("fch1",1400.0,"<d"),
            ("foff",-1.0,"<d"),
            ("nchans",4,"<i"),
            ("nbits",8,"<i"),
            ("nifs",1,"<i"),
        ] + ([("nsamples",samples,"<i")] if with_nsamples else []):
            field(f,key)
            if fmt is None:
                field(f,value)
            else:
                f.write(struct.pack(fmt,value))
        field(f,"HEADER_END")
        f.write(bytes((j*19+11)%256 for j in range(samples*4)))
    return samples


class TestExtractor(unittest.TestCase):
    def test_optional_nsamples(self):
        for with_nsamples in (False,True):
            with self.subTest(with_nsamples=with_nsamples), tempfile.TemporaryDirectory() as d:
                d=Path(d)
                src=d/"source.fil"
                count=make_fil(src,with_nsamples)
                before=src.read_bytes()
                out=d/"windows"
                proc=subprocess.run([sys.executable,str(EXTRACTOR),"--input",str(src),
                                     "--output-dir",str(out),"--windows","4"],
                                    capture_output=True,text=True)
                self.assertEqual(proc.returncode,0,proc.stdout+"\n"+proc.stderr)
                self.assertEqual(src.read_bytes(),before)
                manifest=json.loads((out/"windows_manifest.json").read_text())
                self.assertEqual(sum(w["nsamples"] for w in manifest["segments"]),count)
                self.assertEqual(len(manifest["segments"]),4)
                self.assertTrue(all(Path(w["path"]).exists() for w in manifest["segments"]))
                self.assertFalse(any(out.glob("*.part")))
                self.assertTrue(all(w["expected_bytes"]==Path(w["path"]).stat().st_size
                                    for w in manifest["segments"]))

    def test_does_not_overwrite_existing_results(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d)
            src=d/"source.fil"
            make_fil(src,False)
            out=d/"windows"
            cmd=[sys.executable,str(EXTRACTOR),"--input",str(src),
                 "--output-dir",str(out),"--windows","4"]
            self.assertEqual(subprocess.run(cmd,capture_output=True).returncode,0)
            self.assertNotEqual(subprocess.run(cmd,capture_output=True).returncode,0)


if __name__=="__main__":
    unittest.main()
