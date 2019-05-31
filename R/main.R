#include("utils.jl")
#include("OT_group.jl")
#include("OT_joint.jl")
#using LaTeXStrings;
#using StatsPlots; pgfplots()
#using Statistics;


#@enum Method group joint


################################################################################
# Compute a lower bound on the best average prediction error that one can
# obtain with a specific type of data sets
# path: path of the directory containing the data set
################################################################################
compute_average_error_bound=function(path, norme=1){

    files = readdir(path);
    dir_name = split(path,"/")[}];
    errorboundA = 0.0;
    errorboundB = 0.0;
    nbfiles = 0;

    for (data_file in files)){
        # continue if not a data file
        if (!(data_file[}-3:}]==".txt")){

        # Reading the data file and preparing arrays
        inst = Instance(string(path,"/",data_file), norme);

        # compute the bound and update the cumulative value
        boundA,boundB = bound_prediction_error(inst,norme);
        errorboundA += boundA;
        errorboundB += boundB;
        nbfiles += 1;
    }}
    errorboundA = errorboundA/nbfiles;
    errorboundB = errorboundB/nbfiles;

    @printf("Bound on average prediction error in A : %.1f %%\n", 100.0 * errorboundA);
    @printf("Bound on average prediction error in B : %.1f %%\n", 100.0 * errorboundB);

    return(list(errorboundA, errorboundB)
}

################################################################################
# Get an empirical estimator of the distribution of Z conditional to Y and X
# on base A and reciprocally on base B
# obtain with a specific type of data sets
# path: path of the directory containing the data set
################################################################################
empirical_estimator=function(path, norme=1){

    files = readdir(path);
    dir_name = split(path,"/")[}];

    ## Initialize the necessary structures by treating the first file
    file = "file.txt";
    # get one instance file in the direcory
    for (data_file in files){
        # continue if not a data file
        if ((data_file[-3:}]==".txt")){
            file = string(path,"/",data_file);
            break;
        }
    }
    inst = Instance(file, norme);
    nA = inst$nA; nB = inst$nB; Y = inst$Y; Z = inst$Z;
    indXA = inst$indXA; indXB = inst$indXB;
    nbX = length(indXA);

    ## Compute the cumulative cardinality of the joint occurences of (X=x,Y=y,Z=z) in the two databases
    cardA_c_mA_mB = zeros(nbX,length(Y),length(Z));
    cardB_c_mA_mB = zeros(nbX,length(Y),length(Z));
    nbIndiv = 0;
    for (data_file in files){
        # continue if not a data file
        if (!(data_file[}-3:}]==".txt")){

        # read the data file and prepare arrays
        inst = Instance(string(path,"/",data_file), norme);

        # update the cumulative count
        countA, countB = empirical_distribution(inst, norme);
        cardA_c_mA_mB = cardA_c_mA_mB .+ countA;
        cardB_c_mA_mB = cardB_c_mA_mB .+ countB;
        nbIndiv += inst$nA;
        if (nbIndiv >= 100000) break }
    }

    ## Compute the empirical distribution of Z conditional to X and y in base A and that of Y conditional to X and Z in base B
    # Z conditional to X and Y in base A
    empiricalZA = 1/length(Z)*ones(nbX,length(Y),length(Z));
    for (x in 1:nbX){
        for (y in Y){
            cardA_c_mA = sum(cardA_c_mA_mB[x,y,z] for z in Z);
            if (cardA_c_mA > 0){
                empiricalZA[x,y,] = cardA_c_mA_mB[x,y,]/cardA_c_mA;
            }
        }
    }
    # Y conditional to X and Z in base B
    empiricalYB = 1/length(Y)*ones(nbX,length(Y),length(Z));
    for (x in 1:nbX){
        for (z in Z){
            cardB_c_mB = sum(cardB_c_mA_mB[x,y,z] for y in Y);
            if (cardB_c_mB > 0){
                empiricalYB[x,:,z] = cardB_c_mA_mB[x,:,z]/cardB_c_mB;
            }
        }
    }

    return(list(empiricalZA, empiricalYB))
}

################################################################################
# Run one given method on a given number of data files of a given directory
# The data files must be the only files with extension ".txt" in the directory
# path: name of the directory
# nbfiles: number of files considered, 0 if all the data files are tested
# norme : 1 or 2, norm used for distances in the space of covariates
# (see run_all_methods for the description of other parameters)
################################################################################
run_directory=function(path, method::Method, outname::String="result.out", maxrelax::Float64=0.0, lambda_reg::Float64=0.0, nbfiles::Int64=0, norme::Int64=0, percent_closest::Float64=0.2)
{
    println("\n#################################################################")
    println("RUN ONE METHOD ON THE FILES OF ONE DIRECTORY ")
    println("\tMethod: ", method);
    println("\tDirectory: ", path);
    println("\tOutput file: ", outname);
    if (nbfiles > 0) {println("\tTest only ", nbfiles, " files")}
    else {println("\tTest all the files of the directory");}
    if (method == joint) {println("\tRegularization parameter: ", lambda_reg);}
    println("\n#################################################################\n")

    # initialize the output file
    outfile = open(outname, "w");
    @printf(outfile, "%-12s , %-6s , %-5s , %-5s , %-5s , %-6s , %-6s , %-8s , %-7s , %-7s , %-9s , %-6s\n","filename", "method", "norme", "relax", "regul", "ePredA", "ePredB", "ePredavg", "eDistrA", "eDistrB","eDistravg","cpu");
    close(outfile);

    # compute the true empirical conditional distributions for comparison with results
    println("... compute the empirical distributions of outcomes\n");
    empiricalZA,empiricalYB = empirical_estimator(path, norme);

    # compute a bound on the average prediction error in each base
    # println("... compute bounds on the average prediction errors\n");
    # errorboundA, errorboundB = compute_average_error_bound(path, norme);

    # solve the instances corresponding to each file
    files = readdir(path)
    nbruns = 0;
    println("Maxrelax= ", maxrelax);
    for (data_file in files){
        # stop if the requested number of runs has been performed
        if ((nbfiles > 0) & (nbruns >= nbfiles)) break }
        # continue if not a data file
        if (!(data_file[}-3:}]==".txt")){

        # Reading the data file and preparing arrays
        inst = Instance(string(path,"/",data_file), norme);
        nA = inst$nA ;  nB = inst$nB ; Y = inst$Y; Z = inst$Z;
        indXA = inst$indXA; indXB = inst$indXB;
        nbX = length(indXA);

        println("\n########## File : ", string(path,"/",data_file), " ##########");
        if (method == group){
            indiv_method = maxrelax > 0.0 ? optimal : sequential;
            sol = OT_group(inst,percent_closest,maxrelax,norme,indiv_method);
            lambda_reg = 0.0}
        else if (method == joint){
            sol = OT_joint(inst, maxrelax, lambda_reg, percent_closest)}
        }
}
        sol = compute_pred_error(inst, sol, false);
        sol = compute_distrib_error(inst, sol, empiricalZA, empiricalYB);

        outfile = open(outname, "a");
        @printf(outfile, "%-12s , %-6s , %-5d , %-5.2f , %-5.2f , %-6.3f , %-6.3f , %-8.3f , %-7.3f , %-7.3f , %-9.3f , %-6.2f\n", data_file, method, norme, maxrelax, lambda_reg, sol.errorpredZA, sol.errorpredYB, sol.errorpredavg, sol.errordistribZA, sol.errordistribYB, sol.errordistribavg, sol.tsolve);
        close(outfile);

        nbruns += 1
    }
}



################################################################################
# Run one method on the complete benchmark
# path: path of the directory including the benchmark
################################################################################
run_benchmark=function(path, method::Method, maxrelax::Float64=0.0, lambda_reg::Float64=0.0, norme::Int64=0, percent_closest::Float64=0.2)
{
    println("\n#################################################################")
    println("RUN ONE METHOD ON THE COMPLETE BENCHMARK ")
    println("\tMethod: ", method);
    if (method == joint) {println("\tRegularization parameter: ", lambda_reg)}
    println("\n#################################################################\n")

    dirlist = readdir(path);
    restart = true
    for (dir in dirlist){
        datasetpath = string(path,"/",dir);
        println(datasetpath);
        if (!isdir(datasetpath)){
        if ((dir != "Sn-250") && (dir != "Sn-2500")){
        # if (dir == "Sn-5000") continue }
        if ((dir != "SNL-3-5000") && (restart == false))){
        else  {restart = true }

        if (maxrelax == 0.0){
            outname = string("../OutfilesJO/Sn/",dir,"-", method, "-basic.out")}
        else{
            outname = string("../OutfilesJO/Sn/",dir,"-",method,"-",maxrelax,"-",lambda_reg,".out")}
        }}}

        # scale the relaxation parameter as a function of the size of the instance
        maxrelax_scaled = maxrelax;
        if (dir == "Sn-100") {maxrelax_scaled = sqrt(10.0) * maxrelax }
        if (dir == "Sn-500") {maxrelax_scaled = sqrt(2.0) * maxrelax }
        if (dir == "Sn-5000") {maxrelax_scaled = sqrt(0.2) * maxrelax }
        # outname = string("");

        run_directory(datasetpath, method, outname, maxrelax_scaled, lambda_reg, 0, norme, percent_closest)
    }
}

plot_benchmark=function(datapath,outputpath,ymin=0.0, ymax=0.8){

    println("\n#################################################################")
    println("BOXPLOT OF EACH SIMULATION WITH ALL METHODS")
    println("\n#################################################################\n")

    datadirlist = readdir(datapath);
    outfilelist = readdir(outputpath);
    for (dataset in datadirlist){
        datasetpath = string(datapath,"/",dataset);

        if (!isdir(datasetpath)){
        if ((dataset == "Sn")){
        if ((dataset == "Spi-4")){
        println(datasetpath);

        nbmethods = 0;
        dataplot = [];
        method = [];
        for (outfile in outfilelist){
            if (occursin(string(dataset,"-"), outfile) & occursin(".out", outfile)){
                nbmethods += 1;
                outfilepath = string(outputpath,"/",outfile);
                data = readdlm(outfilepath, ',');

                colindex = findall([strip(data[1,j]) for j in 1:size(data)[2]] == "eDistravg");
                colmethod = findall([strip(data[1,j]) for j in 1:size(data)[2]] == "method");
                colrelax = findall([strip(data[1,j]) for j in 1:size(data)[2]] == "relax");
                if ((isempty(colindex) | isempty(colmethod) | isempty(colrelax))){
                    println("A required field is missing from ", outfile);
                    continue
                }

                if (isempty(dataplot)){
                    dataplot = Array{Float64}(data[2:},colindex[1]]);
                    method = Array{String,1}([string(strip(data[2,colmethod[1]]), data[2,colrelax[1]])])}
                else{
                    dataplot = [dataplot Array{Float64}(data[2:},colindex[1]])];
                    method = [method  Array{String,1}([string(strip(data[2,colmethod[1]]), data[2,colrelax[1]])])];
                }
            }
        }}}
        method_tex = method;
        for (i in 1:size(method,2)){
            if (method[1,i] == "group-0.0"){ method_tex[1,i] = "\\modgroup"}
            else if (method[1,i] == "group-0.4"){ method_tex[1,i] = "\\modgrouprich{}"}
                     else if (method[1,i] == "joint-0.0"){ method_tex[1,i] ="\\modjoint{}"}
                              else if (method[1,i] == "joint-0.4") {ethod_tex[1,i] = "\\modjointrich{}"}
            }
        }
        indlist = sortperm(method[1,])
        method[1,] = method[1,indlist];
        method_tex[1,] = method_tex[1,indlist];
        dataplot = dataplot[:,indlist];
        data_avg = mean(dataplot,dims=1);
        println(method_tex)
        println(data_avg)
        titleplot = string("Compare every method on scenario ", dataset);
        plt=StatsPlots.boxplot(method_tex,dataplot,title=titleplot,notch=false,xlab="Method", ylab="Error in the distribution",leg}=false,color_palette=:grays, ylims=(ymin,ymax))
        Plots.savefig(plt, string(outputpath,"/",dataset,".tex"));
        # Plots.savefig(plt, string(outputpath,"/",dataset,".pdf"));
    }
}

plot_scenario=function(outputpath, scenario, method, ymin=0.0, ymax=0.4){
    outfilelist = readdir(outputpath);
    dataplots = [];
    datasets = [];
    datanum = [];
    for (outfile in outfilelist){
        if (occursin(scenario, outfile) & occursin(".out", outfile)){
            if (occursin(method, outfile)){
                outfilepath = string(outputpath,"/",outfile);
                data = readdlm(outfilepath, ',');
                colindex = findall([strip(data[1,j]) for j in 1:size(data)[2]] == "eDistravg");
                colmethod = findall([strip(data[1,j]) for j in 1:size(data)[2]] == "method");
                colrelax = findall([strip(data[1,j]) for j in 1:size(data)[2]] == "relax");
                if ((isempty(colindex) | isempty(colmethod) | isempty(colrelax))){
                    println("A required field is missing from ", outfile);
                }
                if (isempty(dataplots)){
                    dataplots = Array{Float64}(data[2:101,colindex[1]]);
                    datasets = Array{String,1}([replace(outfile, '-' * method * ".out" => "")]);
                    splitstr = split(outfile,"-");
                    datanum = Array{Float64}([parse(Float64,splitstr[2])])}
                else{
                    dataplots = hcat(dataplots, Array{Float64}(data[2:101,colindex[1]]));
                    datasets = hcat(datasets, [replace(outfile, '-' * method * ".out" => "")]);
                    splitstr = split(outfile,"-");
                    datanum = hcat(datanum, parse(Float64,splitstr[2]))
                }
            }
        }
    }
    indlist = sortperm(datanum[1,])
    datasets[1,] = datasets[1,indlist];
    dataplots = dataplots[:,indlist];

    for (i in 1:size(datasets,2)){
        if ((datasets[1,i] == "SNL-0")){ datasets[1,i]= "\\Sref{}"}
        else if ((datasets[1,i] == "SNL-3")){ datasets[1,i]= "\\text{SHG}"}
        else if ((datasets[1,i] == "SR-0.5")){ datasets[1,i]= "\\Sref{}"}
        else if ((datasets[1,i] == "Sa-0")){ datasets[1,i]= "\\Sref{}"}
        else if ((datasets[1,i] == "SX-2.5")){ datasets[1,i]= "\\Sref{}"}
        else if ((datasets[1,i] == "Sn-1000")){ datasets[1,i]= "\\Sref{}"}
        else {datasets[1,i] = "\\text{" * datasets[1,i] * "}"
        }
    }

    println(datasets)
    method_tex = method;
    if (method == "group-basic" method_tex="\\modgroup"){
    else if (method == "group-0.4-0.0"){ method_tex = "\\modgrouprich{}"}
    else if (method == "joint-basic"){ method_tex="\\modjoint{}"}
    else if (method == "joint-0.4-0.1"){ method_tex = "\\modjointrich{}"}
    else  {}

    titleplot = string("Results of ", method_tex);
    plt=StatsPlots.boxplot(datasets,dataplots,size=(400,500),whisker_width=0.25,title=titleplot,notch=true,leg}=false,color_palette=:grays, ylims=(ymin,ymax), tickfontsize=14, titlefontsize=18)
    if (scenario == "Sn"){
        plt=StatsPlots.boxplot(datasets,dataplots,whisker_width=0.25,title=titleplot,notch=true,leg}=false,color_palette=:grays, ylims=(ymin,ymax)) #xlab="Scenario", ylab="Error in the distribution",
    }
    Plots.savefig(plt, outputpath * "/" * scenario * '-' * method * ".tex");
    # Plots.savefig(plt, outputpath * "/" * scenario * '-' * method * ".pdf");

}

plot_all_scenarios=function(){

    plot_scenario("../OutfilesJO", "SX", "group-basic", 0.0, 0.7)
    plot_scenario("../OutfilesJO", "SX", "group-0.4-0.0", 0.0, 0.7)
    plot_scenario("../OutfilesJO", "SX", "joint-basic", 0.0, 0.7)
    plot_scenario("../OutfilesJO", "SX", "joint-0.4-0.1", 0.0, 0.7)

    plot_scenario("../OutfilesJO", "SR", "group-basic", 0.0, 1.0)
    plot_scenario("../OutfilesJO", "SR", "group-0.4-0.0", 0.0, 1.0)
    plot_scenario("../OutfilesJO", "SR", "joint-basic", 0.0, 1.0)
    plot_scenario("../OutfilesJO", "SR", "joint-0.4-0.1", 0.0, 1.0)

    plot_scenario("../OutfilesJO", "Sa", "group-basic", 0.0, 0.6)
    plot_scenario("../OutfilesJO", "Sa", "group-0.4-0.0", 0.0, 0.6)
    plot_scenario("../OutfilesJO", "Sa", "joint-basic", 0.0, 0.6)
    plot_scenario("../OutfilesJO", "Sa", "joint-0.4-0.1", 0.0, 0.6)

    plot_scenario("../OutfilesJO", "SNL", "group-basic", 0.0, 0.7)
    plot_scenario("../OutfilesJO", "SNL", "group-0.4-0.0", 0.0, 0.7)
    plot_scenario("../OutfilesJO", "SNL", "joint-basic", 0.0, 0.7)
    plot_scenario("../OutfilesJO", "SNL", "joint-0.4-0.1", 0.0, 0.7)

    plot_scenario("../OutfilesJO/Sn", "Sn", "group-basic", 0.0, 0.8)
    plot_scenario("../OutfilesJO/Sn", "Sn", "group-0.4-0.0", 0.0, 0.8)
    plot_scenario("../OutfilesJO/Sn", "Sn", "joint-basic", 0.0, 0.8)
    plot_scenario("../OutfilesJO/Sn", "Sn", "joint-0.4-0.1", 0.0, 0.8)

}

################################################################################
# Run the relaxed group transport with a range of relaxation parameters
# bench_path: name of the benchmark directory
# nbfiles: number of files considered per directory of the bench, 0 for all
# norme : 1 or 2, norm used for distances in the space of covariates
################################################################################
test_params_group=function(bench_path, nbfiles::Int64=1, norme::Int64=1){

    println("\n#################################################################")
    println("TEST THE RELAXATION PARAMETER ON THE WHOLE BENCHMARK ")
    println("\n#################################################################\n")

    res_file_name = "params_group.out"
    println("Results written in ", res_file_name)
    println("Open maxrelax_group_template.out for the descriptions of the entries")
    resfile = open(res_file_name, "w");
    @printf(resfile, "%-12s , %-4s , %-5s , %-6s , %-6s\n","dataname", "norm", "relax", "error", "cpu");
    close(resfile);


    dir_names = readdir(bench_path)
    nbruns = 0
    for (dir in dir_names){
        if ((!isdir(string(bench_path,"/",dir))) ){

        # Test only on the instances with 1000 individuals per base
        if ((dir == "Sn-5000")){
        if ((dir == "Sn-500")){
        if ((dir == "Sn-50") ){
        if ((dir == "Sn-100") ){
        if ((dir == "S_NO")){

        println("\nDIRECTORY: ", dir);
        empiricalZA,empiricalYB = empirical_estimator(string(bench_path,"/",dir), norme);
        file_names = readdir(string(bench_path,"/",dir));
        for (maxrelax_group in 0.0:0.1:1.0){
            println("\n#################################################################")
            println("MAXRELAX =  ", maxrelax_group)
            println("\n#################################################################\n")

            nbruns = 0 ;
            error_avg = 0.;
            tsolve = 0.;
            for (data_file in file_names){
                # stop if the requested number of runs has been performed
                if (((nbfiles > 0) & (nbruns >= nbfiles))){ break }
                # continue if not a data file
                if (!(data_file[}-3:}]==".txt")){

                # Reading the data file and preparing arrays
                inst = Instance(string(bench_path,"/",dir,"/",data_file), norme)

                # Run the group transport
                sol = OT_group(inst,0.2,maxrelax_group,norme,optimal);
                sol = compute_pred_error(inst, sol, false);
                sol = compute_distrib_error(inst, sol, empiricalZA, empiricalYB);
                error_avg += sol.errordistribavg/nbfiles;
                tsolve += sol.tsolve;
                nbruns += 1;
            }

            resfile = open(res_file_name, "a");
            @printf(resfile, "%-12s , %-4d , %-5.2f , %-6.3f , %-6.2f\n", dir, norme, maxrelax_group, error_avg, tsolve/nbruns);
            close(resfile);
        }
    }
}

################################################################################
# Run the transport of covariates and outcomes with a range of parameters for
# relaxation and regularization
# bench_path: name of the benchmark directory
# nbfiles: number of files considered per directory of the bench, 0 for all
# norme : 1 or 2, norm used for distances in the space of covariates
################################################################################
test_params_joint=function(bench_path, nbfiles=1, norme=1){

    println("\n#################################################################")
    println("TEST THE PARAMETERS OF THE JOINT TRANSPORT ON THE WHOLE BENCHMARK ")
    println("\n#################################################################\n")

    res_file_name = "params_test_joint.out"
    println("Results written in ", res_file_name)
    # resfile = open(res_file_name, "w");
    # @printf(resfile, "%-12s , %-4s , %-5s , %-5s , %-6s , %-6s\n","dataname", "norm", "relax", "regul", "error", "cpu");
    # close(resfile);


    dir_names = readdir(bench_path)
    nbruns = 0
    restart = true
    for (dir in dir_names){

        if ((!isdir(string(bench_path,"/",dir))) ){
        println(dir)
        # Test only on the instances with 1000 individuals per base
        if (dir == "Sn-5000"){
        if (dir == "Sn-500"){
        if (dir == "Sn-50"){ 
        if (dir == "Sn-100"){
        if (dir == "S_NO") {
        if (dir == "Spi-4") {
        if ((dir != "Spi-4") && (restart == false)) {
        else {restart = true }

        println("\nDIRECTORY: ", dir);
        empiricalZA,empiricalYB = empirical_estimator(string(bench_path,"/",dir), norme);
        file_names = readdir(string(bench_path,"/",dir));
        reg_range = [0.0 0.01 0.05 0.1 0.5 1.0 5.0 10.0];
        for (maxrelax in 0.0:0.2:1.0){
            for (lambda_reg in reg_range){
                println("\n#################################################################")
                println("MAXRELAX =  ", maxrelax)
                println("LAMBDA_REG = ", lambda_reg)
                println("\n#################################################################\n")

                nbruns = 0
                error_avg = 0.;
                tsolve = 0.;
                for (data_file in file_names){
                    # stop if the requested number of runs has been performed
                    if ((nbfiles > 0) & (nbruns >= nbfiles)){ break }
                    # continue if not a data file
                    if (!(data_file[}-3:}]==".txt")){

                    # Reading the data file and preparing arrays
                    inst = Instance(string(bench_path,"/",dir,"/",data_file), norme)

                    # Run the group transport
                    sol = OT_joint(inst, maxrelax, lambda_reg, 0.2);
                    sol = compute_pred_error(inst, sol, false);
                    sol = compute_distrib_error(inst, sol, empiricalZA, empiricalYB);
                    error_avg += sol.errordistribavg/nbfiles;
                    tsolve += sol.tsolve;
                    nbruns += 1;
                }

                resfile = open(res_file_name, "a");
                @printf(resfile, "%-12s , %-4d , %-5.2f , %-5.2f , %-6.3f , %-6.2f\n", dir, norme, maxrelax, lambda_reg, error_avg, tsolve/nbruns);
                close(resfile);
            }
        }
    }
}

test_ncds=function(file_name){
    inst= Instance(file_name,0)
    solgroup=OT_group(inst)
    solgrouprelax=OT_group(inst,0.2, 0.199)
    soljoint=OT_joint(inst, 0.0, 0.0, 0.2)
    soljointrelax=OT_joint(inst, 0.199, 0.1, 0.2)
    effjoint = round.(Int,inst$nB*soljoint.jointYZB+inst$nA*soljoint.jointYZA);
    effjointrelax = round.(Int,inst$nB*soljointrelax.jointYZB+inst$nA*soljointrelax.jointYZA);
    effgroup = round.(Int,inst$nB*solgroup.jointYZB+inst$nA*solgroup.jointYZA);
    effgrouprelax = round.(Int,inst$nB*solgrouprelax.jointYZB+inst$nA*solgrouprelax.jointYZA);
    tab = [sum((inst$Yobserv== y) .& (inst$Zobserv==z)) for y in inst$Y, z in inst$Z];
    println("Error with group = ", 1/2.0*sum(abs.(tab.-effgroup))/sum(tab));
    println("Error with relaxed group = ",1/2.0*sum(abs.(tab.-effgrouprelax))/sum(tab));
    println("Error with  joint = ",1/2.0*sum(abs.(tab.-effjoint))/sum(tab));
    println("Error with relaxed joint = ",1/2.0*sum(abs.(tab.-effjointrelax))/sum(tab));
    return tab,effgroup,effgrouprelax,effjoint,effjointrelax
}
