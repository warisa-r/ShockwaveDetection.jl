using ShockwaveDetection
using ShockwaveProperties
using Euler2D:Euler2D

flow_data = FlowData("examples/data/obstacle/funky_square.celltape")
point_detect_algo = ImageProcessingShockDetectionAlgo(0.02, :prewitt)
dbscan_algo = DBSCANAlgo(0.02, 10, 10)

# What seems to be happening in the first frame is that
# The shocks seem to be detected at both edges of the square and because the tail of the shock that expands out
# have pretty weak gradients, the algorithm naturally fails to detect the curve above and below the obstacle and fits a line instead.
# Even with finetuning DBSCAN parameters to seperate these two different shock patterns, the
# algorithm that cannot find the tail of the shock will not be able to detect the full curve expanding out of the obstacle
# Decreasing the gradient threshold will increase the shock detects at the edges of the obstacle, making it harder to seperate
# both moving shocks and a line will still be fitted, despite the algorithm having detected more "tail" of the shock

detection = detect(flow_data, point_detect_algo, dbscan_algo)

#plot_shock_fits_over_time(flow_data, detection, false)
create_heatmap_evo_with_shock(flow_data, detection, :density_field, true, false)
