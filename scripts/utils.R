group_colors = c(
    "M"="dodgerblue4",
    "C"="gold",
    "K"="orange",
    "VH"="red4",
    "VL"="red2")
time_colors = grDevices::rainbow(14)
names(time_colors) = c(0, 3, 6, 9, 12, 18, 24, 30, 36, 48, 60, 72, 120, 168)

ann_colors = list(
    Time=time_colors,
    Group=group_colors
    )

replicate_markers = c(15, 17, 19)
names(replicate_markers) = c(1, 2, 3)
ann_markers = list(
    Replicate=replicate_markers)
