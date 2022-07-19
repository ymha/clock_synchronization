This directory includes data recorded by the system described in Section VI of the paper.

  - k            : the rightmost subscript k denotes the repetition number of PTP packet exchange
  
  - offset_k     : clock offset between master and slave nodes (unit: nano-second)
                   offset_k = t_{1,k} - C(t_{1,k})

  - skew_k       : clock skew between master and slave nodes (no unit)

  - timestamps_k : PTP timestamps [t{1,k}, C(t_{2,k}), C(t_{3,k}), t_{4,k}]
                   t_{1,k} is the transmitted time of the k-th PTP sync message recorded by the master node
                   C(t_{2,k}) is the received time of the k-th PTP sync message recorded by the slave node
                   C(t_{3,k}) is the transmitted time of the k-th PTP delay-req message recorded by the slave node
                   t_{4,k} is the received time of the k-th PTP delay-req message recorded by the master node