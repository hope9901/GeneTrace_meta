"""Download public reference inputs atomically; preserve existing files."""
import gzip
import os
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from genetrace.io import sha256, write_json

ROOT = Path(__file__).resolve().parents[1]
BASE = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/'
PINNED = {
    'genomic.fna.gz': '1ec41f951b35f2752d5065edeafe6972455609d4383fb36dae12519a60b27a73',
    'genomic.gff.gz': 'bbcb56dec2d2a00101fb85be4752fbda9895e27023e0c5e08340fc30d0ad5510',
    'protein.faa.gz': '43b97707ae1a08a4b30c200aba6e5c39be47c674321b096478065a2400b1494e',
}


def main():
    output = ROOT / 'data/source'
    output.mkdir(parents=True, exist_ok=True)
    records = []
    for filename in ('genomic.fna.gz', 'genomic.gff.gz', 'protein.faa.gz'):
        destination = output / filename
        url = BASE + 'GCF_000146045.2_R64_' + filename
        if not destination.exists():
            partial = destination.with_suffix('.partial')
            with urllib.request.urlopen(url, timeout=60) as response, partial.open('wb') as handle:
                while True:
                    block = response.read(1024 * 1024)
                    if not block:
                        break
                    handle.write(block)
            with gzip.open(partial, 'rb') as handle:
                while handle.read(1024 * 1024):
                    pass
            os.replace(partial, destination)
        if sha256(destination) != PINNED[filename]:
            raise ValueError('Reference snapshot checksum changed: ' + filename +
                             '. Inspect the upstream update before changing the pinned dataset.')
        records.append(dict(file=filename, url=url, sha256=sha256(destination)))
        print(filename + ': ' + sha256(destination))
    write_json(output / 'sources.json', records)


if __name__ == '__main__':
    main()
