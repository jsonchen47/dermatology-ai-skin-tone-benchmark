# Model Prompt Record

## Claude

The following prompt was recorded for the Claude evaluation, with the input filename adapted to the relevant batch:

> Can you go through this csv sheet called bcc_scc_analysis_claude.csv.
>
> For each of the images listed, use the md5hash column to search for the corresponding image in the finalfitz17k folder (named as `<md5hash>.jpg`). Using Sonnet 4.6, classify what the top 3 skin lesions you think is most likely. List these under the columns cl_dx_1, cl_dx_2, and cl_dx_3. List the percent likelihood that you think it's each of these diagnoses under columns cl_%_1, cl_%_2, and cl_%_3. List any issues you have under cl_notes. Also list the overall likelihood that the lesion is dangerous under column cl_score.

## ChatGPT

The ChatGPT prompt followed the same requested output structure with ChatGPT-specific column prefixes. The exact final prompt text should be inserted here before public release if it can be recovered from the original task record. Do not describe the prompt as verbatim in the manuscript until this is verified.

## Smartphone Applications

Smartphone outputs were collected through each consumer-facing interface. The manuscript Methods should report the tested app version, operating system/device, any required body-site input, and the number of images requiring screenshot re-upload.
