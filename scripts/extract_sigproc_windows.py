#!/usr/bin/env python3
"""Losslessly split a SIGPROC .fil into independently usable sub-filterbanks.

Only tstart changes; nsamples changes when the input header includes it. Original is read-only. Does not overwrite
existing outputs. Writes *.part and only renames on successful verification.
"""
import argparse
import json
import math
import os
from pathlib import Path
import struct

INTS = {'machine_id','telescope_id','data_type','nchans','nbits','nifs',
        'scan_number','barycentric','pulsarcentric','nbeams','ibeam','nsamples'}
DOUBLES = {'tstart','tsamp','fch1','foff','refdm','az_start','za_start',
           'src_raj','src_dej'}
STRINGS = {'source_name','rawdatafile'}

def exact(f, n):
    b=f.read(n)
    if len(b)!=n:
        raise ValueError('Truncated SIGPROC file/header')
    return b

def read_string(f):
    n=struct.unpack('<i',exact(f,4))[0]
    if not (0<=n<=80):
        raise ValueError(f'Invalid SIGPROC keyword/string length: {n}')
    return exact(f,n).decode('utf-8')

def header(f):
    if read_string(f)!='HEADER_START':
        raise ValueError('Not a SIGPROC filterbank')
    fields={}
    locations={}
    while True:
        key=read_string(f)
        if key=='HEADER_END':
            break
        pos=f.tell()
        if key in INTS:
            val=struct.unpack('<i',exact(f,4))[0]
        elif key in DOUBLES:
            val=struct.unpack('<d',exact(f,8))[0]
        elif key in STRINGS:
            val=read_string(f)
        elif key=='signed':
            val=struct.unpack('<b',exact(f,1))[0]
        else:
            raise ValueError(f'Unknown SIGPROC header field: {key}')
        fields[key]=val
        locations[key]=pos
    offset=f.tell()
    f.seek(0)
    return fields,locations,exact(f,offset),offset

def split(src, folder, divisions=4):
    src=Path(src).expanduser().resolve(strict=True)
    folder=Path(folder).expanduser().absolute()
    if src.parent==folder:
        raise ValueError('Refusing to place output next to input without distinct directory')
    with src.open('rb') as inp:
        fields,loc,head,data_offset=header(inp)
        if not {'tstart','tsamp','nbits','nchans'}<=fields.keys():
            raise ValueError(f'Missing required SIGPROC fields: {sorted({"tstart","tsamp","nbits","nchans"}-fields.keys())}')
        bits_per_sample=fields['nchans']*fields['nbits']*fields.get('nifs',1)
        if bits_per_sample%8:
            raise ValueError('This lossless splitter requires integral bytes per time sample')
        bps=bits_per_sample//8
        if not (bps>0 and math.isfinite(fields['tsamp']) and fields['tsamp']>0):
            raise ValueError('Invalid sample geometry')
        payload=src.stat().st_size-data_offset
        if payload<1 or payload%bps:
            raise ValueError('Invalid or truncated data length')
        n=payload//bps
        if 'nsamples' in fields and fields['nsamples'] not in (0,n):
            raise ValueError(f'Header nsamples={fields["nsamples"]}, payload samples={n}; abort')
        if divisions>n:
            raise ValueError('More windows than input time samples')
        plans=[]
        for k in range(divisions):
            start=k*n//divisions
            stop=(k+1)*n//divisions
            output=folder/f'W{k+1}_M62_background.fil'
            partial=Path(str(output)+'.part')
            if output.exists() or partial.exists():
                raise FileExistsError(f'Refusing to overwrite: {output} or {partial}')
            newstart=fields['tstart']+start*fields['tsamp']/86400.0
            plans.append(dict(window=f'W{k+1}',path=str(output),start_sample=start,
                              nsamples=stop-start,tstart_mjd=newstart,
                              duration_seconds=(stop-start)*fields['tsamp'],
                              expected_bytes=len(head)+(stop-start)*bps))
        folder.mkdir(parents=True,exist_ok=True)
        print(f'INPUT {src}: {n} samples; {fields["tsamp"]:.12g} s/sample; {bps} bytes/sample',flush=True)
        for p in plans:
            output=Path(p['path'])
            partial=Path(str(output)+'.part')
            modified=bytearray(head)
            struct.pack_into('<d',modified,loc['tstart'],p['tstart_mjd'])
            if 'nsamples' in loc:
                struct.pack_into('<i',modified,loc['nsamples'],p['nsamples'])
            inp.seek(data_offset+p['start_sample']*bps)
            to_copy=p['nsamples']*bps
            with partial.open('xb') as out:
                out.write(modified)
                while to_copy:
                    buf=inp.read(min(to_copy,16*1024*1024))
                    if not buf:
                        raise IOError('Unexpected EOF while copying filterbank samples')
                    out.write(buf)
                    to_copy-=len(buf)
                out.flush()
                os.fsync(out.fileno())
            if partial.stat().st_size!=p['expected_bytes']:
                raise IOError('Incorrect output size')
            with partial.open('rb') as test:
                hdr,_,_,off=header(test)
            for key,val in fields.items():
                expected=p['tstart_mjd'] if key=='tstart' else p['nsamples'] if key=='nsamples' else val
                if hdr[key]!=expected:
                    raise IOError(f'Header mismatch in {p["window"]}: {key}')
            for sample_index in [0,p['nsamples']//2,p['nsamples']-1]:
                inp.seek(data_offset+(p['start_sample']+sample_index)*bps)
                with partial.open('rb') as test:
                    test.seek(off+sample_index*bps)
                    if exact(test,bps)!=exact(inp,bps):
                        raise IOError(f'Payload mismatch in {p["window"]}, sample {sample_index}')
            os.replace(partial,output)
            print(f'OK {p["window"]}: {p["nsamples"]} samples; {p["duration_seconds"]:.3f}s; {output}',flush=True)
        manifest=folder/'windows_manifest.json'
        with manifest.open('x') as f:
            json.dump(dict(source=str(src),source_header_bytes=data_offset,
                           source_actual_samples=n,source_header=fields,segments=plans),f,indent=2)
            f.write('\n')
        print(f'MANIFEST {manifest}',flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',required=True)
    parser.add_argument('--output-dir',required=True)
    parser.add_argument('--windows',type=int,default=4)
    args=parser.parse_args()
    if not 1<=args.windows<=20:
        parser.error('--windows must be 1..20')
    split(args.input,args.output_dir,args.windows)
