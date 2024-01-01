# Plot_cytoskeletal_orientation
The cytoskeletal orientational order parameter provides a measure of alignment of stress-fibres in a cell. In Buskermolen & Suresh et al [1], we we seeded myofibroblasts on substrates micropatterned with stripes of fibronectin, and observed the alignment of cells and their actin cytoskeleton (i.e. contact guidance) with varying stripe widths. Immunofluorescence images of myofibroblasts on stripes of 10 different widths are shown below (stripe width provided in top left corner, and associated experimental label in the bottom left corner).
![stripes](https://github.com/hamsinisuresh/Plot_cytoskeletal_orientation/assets/46113011/85a05388-b57b-4760-aee6-732c3ece88ba)

Here, we provide the raw data and scripts used to process actin-stained images to visualize local orientation of cytoskeleton. Briefly, we highlight the actin fibres in the green channel of the immunofluorescence images by converting to grayscale and sharpening, followed by edge detection (by convolution with a Laplacian of Gaussian (LoG) kernel) and finally brightness adjustment to highlight the low intensity fibres (panel (a) in the schematic figure below). Panel (b) shows a myofibroblast on the homogeneous substrate (i.e. substrate completely covered with fibronectin (or) infinite-width substrate) in which the actin stress-fibres are detected from “fibreness” estimation, and coloured by their respective orientations with respect to the stripe direction, and by the rotationally-invariant orientation measure used to compute the cytoskeletal orientation parameter.
![figure_s1-01](https://github.com/hamsinisuresh/Plot_cytoskeletal_orientation/assets/46113011/44951e16-7243-4b26-920c-56766510db74)

Raw and processed images are stored in folders named after the experimental labels of the asociated stripes (for example, images from 50 micron stripe are in folder SLB1). Here are some processed images of myofibroblasts from fibronectin stripes of three different widths, colored by the rotationally-invariant orientation measure.
![figure4ab-01](https://github.com/hamsinisuresh/Plot_cytoskeletal_orientation/assets/46113011/faa509b4-8939-448a-90e4-59aad25de986)


For additional details of image processing and calculation of cytoskeletal orientation parameter, please refer to Supplementary Information 1 in:\
[1] Buskermolen & Suresh et al, **Entropic forces drive cellular contact guidance.** Biophysical Journal, 116, 1994–2008, 2019.

If you find this resource useful, please consider citing us.
