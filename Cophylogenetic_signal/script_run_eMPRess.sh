
#####  SCRIPT FOR RUNNING  EMPRESS ON SIMULATIONS WITH VICARIANCE 


### RUN EMPRESS TO GET THE RECONCILIATIONS

export PYTHONUSERBASE=$WORK/cophylo_signal/emPress/script/.local_py3.7.10
export PATH=$PYTHONUSERBASE/bin:$PATH

module load python/3.7.10

cd $WORK/cophylo_signal/emPress/script/empress-v1.2.1

data_path=$WORK"/cophylo_signal/simulated_data_geo/"

name=1
d_rate=4
t_rate=1
l_rate=1

for seed in {1..250}
do
    if test -f $data_path"/../emPress/results_geo/reconciliation_d"$d_rate"_t"$t_rate"_l"$l_rate"_simul_geo_"$name"_"$seed".csv"; then
        echo $seed
    else
        echo $seed

        python empress_cli.py reconcile $data_path"/tree_hosts_only_int_simul_geo_"$name"_"$seed".tre" $data_path"/tree_parasites_only_int_simul_geo_"$name"_"$seed".tre" $data_path"/links_simul_geo_"$name"_"$seed".txt" -d $d_rate -t $t_rate -l $l_rate --csv $data_path"/../emPress/results_geo/reconciliation_d"$d_rate"_t"$t_rate"_l"$l_rate"_simul_geo_"$name"_"$seed".csv"
    fi
done



### RUN EMPRESS TO GENERATE THE P-VALUES OF THE RECONCILIATION

export PYTHONUSERBASE=$WORK/cophylo_signal/emPress/script/.local_py3.7.10
export PATH=$PYTHONUSERBASE/bin:$PATH

module load python/3.7.10

cd $WORK/cophylo_signal/emPress/script/empress-v1.2.1

data_path=$WORK"/cophylo_signal/simulated_data_geo/"

name=1
d_rate=4
t_rate=1
l_rate=1

for seed in {1..250}
do
    if test -f $data_path"/../emPress/results_geo/res_d"$d_rate"_t"$t_rate"_l"$l_rate"_simul_geo_"$name"_"$seed".svg"; then
        echo $seed
    else
        echo $seed
        python empress_cli.py p-value $data_path"/tree_hosts_only_int_simul_geo_"$name"_"$seed".tre" $data_path"/tree_parasites_only_int_simul_geo_"$name"_"$seed".tre" $data_path"/links_simul_geo_"$name"_"$seed".txt" -d $d_rate -t $t_rate -l $l_rate --n 1000 --outfile $data_path"/../emPress/results_geo/res_d"$d_rate"_t"$t_rate"_l"$l_rate"_simul_geo_"$name"_"$seed".svg"
    fi
done


# Use the same script with other cost values 

# Use the same script for simulations of trait matching and codiversification (only the path and the names of the input files need to be updated)





