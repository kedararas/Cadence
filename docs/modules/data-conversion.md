# Data Conversion

Imports raw high-speed optical-mapping data into the CADENCE `cmos_all_data`
format.

- **Input:** folders of raw recordings — SciMedia `.gsd/.gsh` or TIFF stacks.
- **Output:** one `cmos_all_data.mat` per recording (per-camera image stacks,
  preview images, `acqFreq`, and the analog/pacing channel when present).
- **Next step:** [Signal Conditioning](signal-conditioning.md).

!!! tip "Video"
    TODO: embed the Data Conversion tutorial.

## Workflow
1. **Select the input directory** containing the raw recordings (`.TIFF` / `.GSD`).
   Select one or more folders to process sequentially.
2. **Set the sampling frequency (Hz)** for the recordings.
3. **Select the output directory.**
4. **Execute the conversion.** Each recording is written as `<name>.mat`.

## Notes & options
- TODO: list supported camera/file formats and any naming conventions.
- For SciMedia dual-camera `.gsh`, the analog/pacing channel is read from camera 1.
- TODO: note the default sampling frequency fallback.

## Troubleshooting
- *Empty/garbage thumbnail:* TODO.
- *"No .tif or .gsh files found":* the selected folder has no supported raw files.
