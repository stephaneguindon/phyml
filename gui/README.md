Requirements
---------------------

 * Java (version = 6) 
 * maven (version >= 3)
 
Downloads
----------------

I've set up a job on our CI Jenkins instance, the latest executable jar can be downloaded from:

https://code.ceres.auckland.ac.nz/jenkins/job/phyml-grid-SNAPSHOT/
 
Importing into Eclipse
-------------------------------

 * File -> Import... -> Maven -> Existing Maven Projects
 * Next
 * Browse to this directory
 * Finish
 
 
Building outside of IDE
---------------------------------
 
    mvn clean install
 

Running commandline client
-----------------------------------------

    java -cp target/phyml-grid-binary.jar nz.org.nesi.phyml.Client --help    # print usage
	java -cp target/phyml-grid-binary.jar nz.org.nesi.phyml.Client -f <path_to_phyml_inputfile> -w 10h -p <phyml_cli_parameters> -c <cpus>

	
Running gui
-----------

    java -jar target/phyml-grid-binary.jar
