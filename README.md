# TR-23

# MiniMutes/EH Data Processor

The `EH Data Processor` is a Python script developed to process and analyze data derived from the [Expansion Hunter](https://github.com/Illumina/ExpansionHunter) tool. It is designed to handle ndjson formatted output files and can compute relevant statistical metrics for further analysis.

## Getting Started

These instructions will help you get a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

- Python 3.6+
- Pandas library
- Argparse library
- Scipy library

You can install any missing packages using pip:

```bash
pip install -r requirements.txt
```

### Usage

The script can be run from the command line using the following format:

```bash
python EH_Processor.py --EHD <EHD directory> --manifest <manifest file> [other options]
```

- `--EHD`: (required) Directory containing ndjson files to be processed
- `--manifest`: (required) CSV file listing paired files to be processed. The CSV file should contain the following columns: icgc_donor_id, control_object_id, case_object_id, sex
- `--disease`: (optional) Name of the disease. If not provided, the name of the `--EHD` argument will be used as the disease name
- `--pvals`: (optional) Boolean flag indicating whether to calculate p-values (default: True)

The script will output a tidy DataFrame and optionally calculate p-values if the `--pvals` flag is set to `True`.

### Example

```bash
python EH_Processor.py --EHD /path/to/ndjson_files --manifest manifest.csv --disease "Disease Name" --pvals True
```

## Authors

- Lucas Finot
- 
## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
