MSIAC: A Julia package for SIAC filtering
========
_________________________________________________
	Version: 1.0
	For bugs / issues, email: xulia.docampo@pm.me
_________________________________________________


## Installing Julia (skip if you already have it)

   1. ****Download &  install**** [Julia](https://julialang.org/downloads/)
   2. ****Set up an alias to the Julia folder****:
      - ****Linux / bash****: 
      ```bash 
      vi ~/.bash_aliases 
      alias julia="/etc/julia-1.7.2/bin/julia"  (change /etc/julia-1.7.2/  to your path to Julia folder)
      ```
      - ****iOS / zsh****: 
      ```zsh
      vi ~/.zshrc 
      alias julia="/etc/julia-1.7.2/bin/julia"  (change /etc/julia-1.7.2/  to your path to Julia folder)
      ```
   3. ****Check that it works. Source the file (or open a new terminal) and type $julia on terminal****

## Installing the MSIAC package
   1. ****Set up an environment variable****
	This is optional. Useful to call MSIAC from anywhere in your machine but you can skip this and go to step 2. 
      - ****Linux / bash****
      ```bash
      vi ~/.bashrc
      export MSIAC=$your_root_to_MSIAC_folder
      ```
      - ****iOS / zsh****
      ```zsh
      vi ~/.zshrc
      export MSIAC=$your_root_to_this_folder
      ```
   2. ****Load & test the MSIAC Package****
      ```julia
      (1) open julia session inside your MSIAC folder
         $cd MSIAC
         $julia

         NOTE: if you are in an editor, make sure the REPL is indeed in MSIAC:
         julia> cd(ENV["MSIAC"])  #(ENV[""] looks for environment variables)

      (2) Enter Pkg environment
         julia> ]
            (you should now see something like @(v1.7) pkg>)

      (3) Activate MSIAC package
         (@v1.6) pkg> activate .   # The dot (.) tells Julia to look in local folder
            (you should see that (@v1.7) pkg> became (MSIAC) pkg> )

      (4)	Load package
            (MSIAC) pkg> instantiate
            ** If you get errors, try running
                        pkg> resolve ; pkg> update

      (5) Test the package (Optional)
         (MSIAC) pkg>  test
            (If all went well, the screen should show 4 (6) test PASSED )
      ```

	***NOTE****: WriteVTK and VTKDataIO packages show some warnings during pre-compilation. However, they still work and allow reading .vtu & .vtk files correctly. 

  3. ****Known issue****
	For some configurations with iOS, building the package using the Project.toml may not work. If the installation failed, checkout README_pkg.md file


## Running the MSIAC tool on Distributed Memory
	Julia is more efficient using DistributedArrays than the default parallelization used by MSIAC: Threads.@threads
	
	1. The file MSIAC_HPC.jl assumes ther is an environment variable with the MSIAC path
		@everywhere push!(LOAD_PATH,ENV["MSIAC"])
	
	***NOTE****:  the everywhere command tells julia to include MSIAC module in all workers. Otherwise, it will only be loaded to the main worker ( # 1 )
	
	   If you haven't set it, manually change ENV["MSIAC"] to your path to MSIAC:
		@everywhere push!(LOAD_PATH,absolute_path_to_MSIAC)
	
	
	2. Start a session with several workers (CPUs) 
		$ julia --p number_of_CPUS
	
	3. Include the MSIAC_HPC.jl file at the top of your script 
		using DelimitedFiles, Printf, Distributed
		@everywhere include(ENV["MSIAC"]*"/src/MSIAC_HPC.jl") 	
	
	4. Run filter_data_hpc(args...) 
	5. This [file](https://github.com/XuliaDS/tree/main//MSIAC/tutorials/msiac_hpc.jl) is an example of a filter + distributed 

# A quick guide on using the REPL
   1. ****Running a script****
      ```bash
         open Julia session
         $julia
         julia> include("myfile.jl")  (you should be at the folder where myfile lives)
            Leaving a session open means that successive runs on the script are faster as it does not recompile:

         julia> include("myfile.jl")   #1st time: slower
         julia> include("myfile.jl")   #2nd time: faster
      ```
   2. ****Checkout the [tutorials](https://github.com/XuliaDS/tree/main/MSIAC/tutorials) folder to see examples of how to run the MSIAC tool****
