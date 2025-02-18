# DRAGEN-GATK Workflow Usage Guide

## Configuring the Workflow
This workflow is executed via a user interface (UI). Follow these steps:

### Step 1: Upload Input Files
- **Sample Sheet (CSV)**: Upload a structured CSV file containing sample metadata.
- **Reference Genome (FASTA)**: Provide a `.fa` or `.fasta` file, optionally compressed (`.gz`).
- **Known Sites VCF (optional)**: Upload a VCF file for base quality recalibration.

### Step 2: Configure Settings
Within the UI, set parameters:
- **Input Paths**: Ensure all uploaded files are linked correctly.
- **Reference Genome**: Select the uploaded reference file.
- **Output Directory**: Choose where results will be saved.

### Step 3: Start the Workflow
- Click the **Run** button to initiate processing.
- The UI will display logs and progress updates.

### Step 4: Review Results
After completion:
- **Aligned Reads**: Found in the `aligned_reads/` directory.
- **Recalibrated BAMs**: Located in `recalibrated/` if recalibration was enabled.
- **Variant Calls (VCF)**: Available in `variants/`.

## Troubleshooting
- Ensure files are correctly formatted before uploading.
- If an error occurs, check logs via the UI.

For further assistance, refer to the official documentation.
