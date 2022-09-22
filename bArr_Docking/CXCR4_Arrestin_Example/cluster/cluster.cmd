$ROS/source/bin/cluster.linuxgccrelease \
    -database $R/database/  \
    -in:file:fullatom \
    -cluster:skip_align \
    @cluster_flags \
    -in:file:silent ../out*silent \
    -in:file:silent_struct_type binary \
    -cluster:population_weight 1 \
    -cluster:sort_groups_by_energy \
    -out:file:silent_struct_type binary \
    -out:file:silent out_cl.silent \
    > cl.out \
    2> cl.err
