# pip install numpy
# pip install pandas
# pip install matplotlib
# pip install sklearn
# apt install python-tk
# install.packages("caTools")
# install.packages("rjson")
# [KernelRidge, AdaBoost, BernoulliNB, MLP, RandomForest, LinearSVC, SVC, knn]

flowps_main = function(
  train_file,
  out_dir=NULL,
  clf='LinearSVC', 
  clf_args=list(),
  min_surround=0,
  max_surround=NULL,
  min_neighbours=32,
  max_neighbours=NULL
){

  prediction.R = function(){
    # specify files:
    predictions_file = file.path(out_dir, "predictions.csv")
    output_file = file.path(out_dir, "y_flowps_score.csv")

    # specify confidence parameter for prediction-accountable set S:
    p = 0.90

    data = read.csv(train_file)
    names = as.vector(data[, 1])
    y = as.vector(data[, 2])  # clinical response
    nsamples = length(y)

    predictions = read.table(predictions_file, sep = ",")

    ylim1 = c(min_surround, max_surround)
    xlim1 = c(min_neighbours, max_neighbours)

    xspan = max_neighbours - min_neighbours + 1
    yspan = max_surround - min_surround + 1

    AUC = array(dim = c(xspan, yspan))
    flowps_score = array(0.0, dim = c(nsamples))

    for(isample in 1:nsamples){
      print(paste0("isample: ", isample))
      y0 = y[-isample]

      # calculates AUC(m, k) for all but one (i.e. isample) samples
      for(ycur in ylim1[1]:ylim1[2]){
        for(xcur in xlim1[1]:xlim1[2]){
          row = xspan * (ycur - ylim1[1]) + (xcur - xlim1[1]) + 1
          scores = as.vector(t(predictions[row,]))
          scores0 = scores[-isample]

          AUC0 = caTools::colAUC(scores0, y0)
          if(cor(scores0, y0) < 0){
            AUC0 = 1 - AUC0
          }

          AUC[xcur - xlim1[1] + 1, ycur - ylim1[1] + 1] = AUC0
        }
      }

      # specifies the prediction-accountable set S for sample isample
      good_row_col = which(AUC > p * max(AUC), arr.ind = TRUE)
      ngood = nrow(good_row_col)
      if(ngood == 0){
        stop("Too few samples.")
      }

      # calculates FloWPS predictions
      for(igood in 1:ngood){
        k = good_row_col[igood, 1] + min_neighbours - 1
        m = good_row_col[igood, 2] - min_surround - 1
        row = xspan * (m - ylim1[1]) + (k - xlim1[1]) + 1
        score = as.vector(t(predictions[row,]))
        flowps_score[isample] = flowps_score[isample] + score[isample]
      }
      flowps_score[isample] = flowps_score[isample] / ngood
    }

    y_flowps_score = cbind(names, y, flowps_score)
    write.table(y_flowps_score, output_file, col.names = TRUE, row.names = FALSE, sep = ",")
  }
  if(is.null(max_neighbours)){
    max_neighbours = length(as.vector(read.csv(train_file)[,2]))-1
  }
  if(is.null(max_surround)){
    max_surround = floor((max_neighbours-1)/2)
  }
  if(is.null(out_dir)){
    out_dir = file.path("./results", basename(train_file))
  }
  system2("python", c(
    system.file("flowps.py", package="flowpspkg"),
    "--train-file", train_file,
    "--out-dir", out_dir,
    "--clf", clf,
    "--clf-args", shQuote(rjson::toJSON(clf_args)),
    "--min-surround", min_surround,
    "--max-surround", max_surround,
    "--min-neighbours", min_neighbours,
    "--max-neighbours", max_neighbours
  ))
  prediction.R()
}
