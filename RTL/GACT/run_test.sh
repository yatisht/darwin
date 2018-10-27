set -x
curr_dir=$PWD
mkdir -p results
cd results/ 
for i in 1 2 3 4 5 6 7 8 9 10;
do
    echo '' > test.${i}.out
    echo "sed -i -e 's/ref\.\(.*\)\.hex/ref\.${i}\.hex/g' ${curr_dir}/tb_GACTTop.v" | sh -e
    echo "sed -i -e 's/query\.\(.*\)\.hex/query\.${i}\.hex/g' ${curr_dir}/tb_GACTTop.v" | sh -e 
    echo "sed -i -e 's/test\.\(.*\)\.vcd/test\.${i}\.vcd/g' ${curr_dir}/tb_GACTTop.v" | sh -e
    iverilog ${curr_dir}/tb_GACTTop.v ${curr_dir}/GACTTop.v ${curr_dir}/BRAM.v ${curr_dir}/FIFOWithCount.v ${curr_dir}/SmithWatermanArray.v ${curr_dir}/SmithWatermanPE.v ${curr_dir}/Ascii2Nt.v ${curr_dir}/BTLogic.v ${curr_dir}/DP_BRAM.v ${curr_dir}/Nt2Param.v ${curr_dir}/mux_1OfN.v
    ./a.out >> test.${i}.out
done
cd $curr_dir
echo 'Results in folder '${curr_dir}'/results/'
