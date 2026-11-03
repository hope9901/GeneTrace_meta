import argparse
import json
from .engine import investigate


def main():
    parser = argparse.ArgumentParser(description='GeneTrace local-locus evidence baseline')
    parser.add_argument('--dataset', required=True, help='Prepared directory with input manifest')
    parser.add_argument('--output', required=True, help='New evidence output directory (must not exist)')
    args = parser.parse_args()
    report = investigate(args.dataset, args.output)
    print(json.dumps(dict(verdict=report['verdict'], mode=report['execution_mode'], output=args.output), indent=2))


if __name__ == '__main__':
    main()
