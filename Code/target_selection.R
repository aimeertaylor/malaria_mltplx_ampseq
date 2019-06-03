#=============================================================
# Previous target selection functions 
# panel_strategy: 'Sample at random', 'Maximise effective cardinality', 
# 'Maximise log(distance) * effective cardinality', 'Maximise distance'
# assuming data has effective cardinality 
#=============================================================
target_selection_static_dt = function(data_, m, panel_strategy){
  
  Mmax = nrow(data_)
  
  # Create distances
  data_$dt <- c(diff(data_$pos), Inf)
  pos_change_chrom <- 1 + which(diff(data_$chrom) != 0) # find places where chromosome changes
  data_$dt[pos_change_chrom-1] <- Inf
  
  # Check distance and effective cardinality are the same magnitude
  data_$Med_Keff*log(data_$dt)
  
  # Order by distance, eff. cardinality, log(distances) * eff. cardinality
  Order_dt = sort.int(log(data_$dt), decreasing = T, index.return = T)$ix
  Order_Keff = sort.int(data_$Med_Keff, decreasing = T, index.return = T)$ix
  Order_logdt_Keff = sort.int(data_$Med_Keff*log(data_$dt), decreasing = T, index.return = T)$ix
  
  if(panel_strategy == "Maximise distance"){
    inds = Order_dt[1:m]} 
  if(panel_strategy == "Maximise effective cardinality"){
    inds = Order_Keff[1:m]}
  if(panel_strategy == "Maximise log(distance) * effective cardinality"){
    inds = Order_logdt_Keff[1:m]}
  if(panel_strategy == 'Sample at random'){
    inds = sample(Mmax, m, replace = F)}
  
  return(data_[sort(inds),])
}




#=============================================================
# New target selection function: can we do better? 
# http://www1.cmc.edu/pages/faculty/MHuber/Research/talks/huber_talk_2011h.pdf
#=============================================================
target_selection_dynmic_dt = function(data_, m, R = 1000){
  
  Mmax = nrow(data_)
  
  #========================================================
  # First step inspired by Matérn Type III process
  #========================================================
  # Sample with probability proportional to diversity
  # Each naturally has a birthday in [0, ∞) - its position
  ind_propose = sort(sample(Mmax, m, prob = data_$Med_Keff, replace = F)) 
  
  # Run time forward, only allowing birth if not within R of older 
  ind_acceptd = ind_propose[1] 
  for(t in 2:m){
    ind_b4 = tail(ind_acceptd,1)
    data_t = data_[c(ind_b4,ind_propose[t]),] # subset data
    dt <- diff(data_t$pos) # calculate dt
    dt[diff(data_t$chrom) != 0] <- Inf # Make inf if chrom changed
    if(dt < R){
      next()
    } else {
      ind_acceptd = c(ind_acceptd, ind_propose[t])
    }
  }
  
  #========================================================
  # Second stepp inspired by Strauss process
  #========================================================
  m_remain = m - length(ind_acceptd)
  
  while(m_remain > 0){
    
    ind_remain = (1:nrow(data_))[-ind_acceptd]
    ind_propose = sort(sample(ind_remain, m_remain, prob = data_$Med_Keff[ind_remain], replace = F)) 
    
    for(t in 1:m_remain){
      
      diffs = ind_acceptd - ind_propose[t]
      diff_sign_diffs = diff(sign(diffs))
      
      if(all(diff_sign_diffs == 0)){
        # If the all accepted are larger
        if(all(diffs > 0)){
          # subset data
          data_t = data_[c(ind_propose[t],ind_acceptd[1]),] 
          dt <- diff(data_t$pos) # Two distances 
          dt[diff(data_t$chrom) != 0] <- Inf 
          if(dt < R){ # If close to either of its neighbours, exclude
            next()} else {
              ind_acceptd = c(ind_propose[t], ind_acceptd)}
        } 
        
        # If the all accepted are smaller
        if(all(diffs < 0)){
          # subset data
          data_t = data_[c(ind_propose[t],ind_acceptd[1]),] 
          dt <- diff(data_t$pos) # Two distances 
          dt[diff(data_t$chrom) != 0] <- Inf 
          if(dt < R){ # If close to either of its neighbours, exclude
            next()} else {
              ind_acceptd = c(ind_acceptd, ind_propose[t])}
        } 
      } else {
        # Find closest neighbours in ind_accepted
        n1_ind = which(diff_sign_diffs != 0)
        n2_ind = n1_ind + 1
        
        # subset data
        data_t = data_[c(ind_acceptd[n1_ind],ind_propose[t],ind_acceptd[n2_ind]),] 
        dts <- diff(data_t$pos) # Two distances 
        dts[diff(data_t$chrom) != 0] <- Inf 
        if(any(dts < R)){ # If close to either of its neighbours, exclude
          next()} else {
            ind_acceptd = c(ind_acceptd, ind_propose[t])}
      }
      ind_acceptd = sort(ind_acceptd)
    }
    m_remain = m - length(ind_acceptd)
    print(m_remain)
  }
  return(data_[sort(ind_acceptd),])
}