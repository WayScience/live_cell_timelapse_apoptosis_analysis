# Mean Average Precision (mAP) Analysis for Wave 1
This folder contains two different analysis.
Each analsysis uses mAP as a metric to answer each question.
1. Does the number of channels effect image-based profiles across treatments that induce different mechanisms of action (MoA)?
2. Does the number of feilds of view (FOV) effect image-based profiles across treatments that induce different mechanisms of action (MoA)?

To test this we used the Wave 1 pyroptosis data image-based profiles.

To adequately answer question 1 we run mAP for each channel individually, and with a leave one channel out approach.
This effectively allows us to see if the number of channels effects the image-based profiles across treatments that induce different mechanisms of action (MoA) via any changes in mAP over time.

To adequately answer question 2 we make an assumptions that the number of fields of view (FOV) is a proxy for the number of cells in the image.
We then run mAP 100 times for each subset of a cell population percentage.
We run the mAP on 10%, 20%, 30%, 40%, 50%, 60%, 70%, 80%, 90%, and 100% of the cell population's aggregated image-based profiles.
This allows us to asses the "stability" of the image-based profiles' mAP score given the differences in the number of fields of view (FOV) or cells in the image.
