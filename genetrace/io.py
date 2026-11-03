import gzip
import hashlib
import json
from pathlib import Path
from urllib.parse import unquote


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def read_fasta(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    records = {}
    name, chunks = None, []
    with opener(path, 'rt') as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    records[name] = ''.join(chunks).upper()
                name, chunks = line[1:].split()[0], []
                if name in records:
                    raise ValueError('Duplicate FASTA identifier: ' + name)
            elif name is None:
                raise ValueError('Sequence before FASTA header')
            else:
                chunks.append(line)
    if name is not None:
        if name in records:
            raise ValueError('Duplicate FASTA identifier: ' + name)
        records[name] = ''.join(chunks).upper()
    if not records or any(not seq for seq in records.values()):
        raise ValueError('Empty FASTA record or file')
    return records


def write_fasta(path, records):
    with open(path, 'w') as handle:
        for name, seq in records.items():
            handle.write('>' + name + '\n')
            for start in range(0, len(seq), 80):
                handle.write(seq[start:start + 80] + '\n')


def read_gff(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    features = []
    with opener(path, 'rt') as handle:
        for line in handle:
            if line.startswith('##FASTA'):
                break
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                raise ValueError('GFF3 requires nine columns')
            attrs = dict((key, unquote(value)) for key, value in
                         (part.split('=', 1) for part in fields[8].split(';') if '=' in part))
            start, end = int(fields[3]), int(fields[4])
            if start < 1 or end < start:
                raise ValueError('Invalid GFF3 coordinates')
            features.append(dict(seqid=fields[0], type=fields[2], start=start,
                                 end=end, strand=fields[6], attrs=attrs, fields=fields))
    return features


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')
