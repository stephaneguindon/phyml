package phyml;

import nz.org.nesi.phyml.swing.GridPanel;
import nz.org.nesi.phyml.swing.PhyMLJobSubmit;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;


/**
 * Class callPHYML provides a method to call the external program PhyML. The
 * parameters are passed in from an external source. The appropriate version of
 * PHYML correlating to the OS (Windows, OSX, Linux 32-bit, Linux 64-bit) with
 * flags corresponding to the parameters is called after determination of the
 * user's OS.
 *
 * @author Vicky Fan, Christoph Knapp
 * @version 14-June-2012
 *
 */
public class callPHYML {

	public static final Logger myLogger = LoggerFactory.getLogger(callPHYML.class);

	/**
	 * Method to call the program PHYML via command line. Note: if unsure of
	 * what the parameters exactly are, please read the PHYML documentation. The
	 * parameters are passed from the GUI.
	 * @param nesiSubmit
	 *
	 * @param seqFileName
	 *            filename of nucleotide/aa file (phylip format)
	 * @param dataType
	 *            nucleotide or amino acid nt/aa
	 * @param sequential
	 *            sequential/interleaved format (file)
	 * @param nDataSets
	 *            number of datasets to analyse
	 * @param parsimony
	 *            use parsimony starting tree
	 * @param bootstrap
	 *            number of bootstrap
	 * @param subModel
	 *            substitution model name
	 * @param aaRateFileName
	 *            aa rate filename if under 'custom' model.
	 * @param freq
	 *            n or aa frequencies
	 * @param tstvRatio
	 *            transition/transversion ratio
	 * @param pInvar
	 *            proportion invariable sites
	 * @param nSubCat
	 *            number of substitution rate categories
	 * @param gamma
	 *            value of gamma shape parameter
	 * @param useMedian
	 *            middle of each substitution rate class in the discrete gamma
	 *            distribution
	 * @param freeRates
	 *            alternative to discrete gamma model
	 * @param codePos
	 *            1,2,4 coding position for estimation
	 * @param search
	 *            tree topology search operation option
	 * @param treeFile
	 *            user tree file
	 * @param params
	 *            specific parameter optimisation
	 * @param randStart
	 *            set initial tree to random
	 * @param nRandStart
	 *            number of initial random trees
	 * @param randSeed
	 *            number used to seed random number generator
	 * @param printSitelnl
	 *            print likelihood for each site in file * phyml lk.txt.
	 * @param printTrace
	 *            print trace of phylogeny explored in * phyml trace.txt
	 * @param runID
	 *            append ID at end of each PhyML output
	 * @param noMemoryCheck
	 *            check with user if want to run even if lots of memory required
	 * @param noJcolAlias
	 *            preprocess each alignment
	 * @param contrainedLens
	 *            find the branch multiplier with input tree with branch length
	 *            provided
	 * @param constFile
	 *            file name topological constraints
	 * @param quiet
	 *            runs PhyML in quiet mode
	 *
	 * @return void
	 *
	 */
	public static void callPhyML(GridPanel gridPanel, boolean nesiSubmit, String seqFileName, String dataType,
			String sequential, String nDataSets, String parsimony,
			String bootstrap, String subModel, String aaRateFileName,
			String freq, String tstvRatio, String pInvar, String nSubCat,
			String gamma, String useMedian, String freeRates, String codePos,
			String search, String treeFile, String params, String randStart,
			String nRandStart, String randSeed, String printSitelnl,
			String printTrace, String runID, String noMemoryCheck,
			String noJcolAlias, String contrainedLens, String constFile,
			String quiet) {

		// format the parameters such that PhyML can be executed
		String command = formatCommand(nesiSubmit, seqFileName, dataType, sequential,
				nDataSets, parsimony, bootstrap, subModel, aaRateFileName,
				freq, tstvRatio, pInvar, nSubCat, gamma, useMedian, freeRates,
				codePos, search, treeFile, params, randStart, nRandStart,
				randSeed, printSitelnl, printTrace, runID, noMemoryCheck,
				noJcolAlias, contrainedLens, constFile, quiet);

		System.out.println("command: " + command);

		try {
			if (command.equals("")) {
				StandardOutPanel
						.setInput("ERROR: Your operating system is not supported...");
				return;
			}
			// extract from the command line string the place where to put the
			// executable
			// file/binaries and call "UnpackExecutableFile" to unpack the right
			// fies there.
			//UnpackExecutableFile.start(command.split(" ")[0]);

            StandardOutPanel.clearPanel();

			if(!nesiSubmit){
				CmdExec exec = new CmdExec(command);
				exec.start();
			}else{
				final PhyMLJobSubmit job = new PhyMLJobSubmit(gridPanel.getServiceInterface(), command, seqFileName);

				Thread t = new Thread() {
					public void run() {
						try {
							job.submit();
						} catch (Exception e) {
							e.printStackTrace();
						}
					};
				};

				t.start();


//				SwingClient frame = new SwingClient(command,seqFileName);
//				frame.run();
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	} // end callPhyML

	/**
	 * Method to return a correctly formated string of the PhyML program
	 * (version) path with the parameters formatted with the correct flags
	 * @param nesiSubmit
	 *
	 * @param seqFileName
	 *            filename of nucleotide/aa file (phylip format)
	 * @param dataType
	 *            nt/aa
	 * @param sequential
	 *            sequential/interleaved format (file)
	 * @param nDataSets
	 *            number of datasets to analyse
	 * @param parsimony
	 *            use parsimony starting tree
	 * @param bootstrap
	 *            number of bootstrap
	 * @param subModel
	 *            substitution model name
	 * @param aaRateFileName
	 *            aa rate filename if under 'custom' model.
	 * @param freq
	 *            n or aa frequencies
	 * @param tstvRatio
	 *            transition/transversion ratio
	 * @param pInvar
	 *            proportion invariable sites
	 * @param nSubCat
	 *            number of substitution rate categories
	 * @param gamma
	 *            value of gamma shape parameter
	 * @param useMedian
	 *            middle of each substitution rate class in the discrete gamma
	 *            distribution
	 * @param freeRates
	 *            alternative to discrete gamma model
	 * @param codePos
	 *            1,2,4 coding position for estimation
	 * @param search
	 *            tree topology search operation option
	 * @param treeFile
	 *            user tree file
	 * @param params
	 *            specific parameter optimisation
	 * @param randStart
	 *            set initial tree to random
	 * @param nRandStart
	 *            number of initial random trees
	 * @param randSeed
	 *            number used to seed random number generator
	 * @param printSitelnl
	 *            print likelihood for each site in file * phyml lk.txt.
	 * @param printTrace
	 *            print trace of phylogeny explored in * phyml trace.txt
	 * @param runID
	 *            append ID at end of each PhyML output
	 * @param noMemoryCheck
	 *            check with user if want to run even if lots of memory required
	 * @param noJcolAlias
	 *            preprocess each alignment
	 * @param contrainedLens
	 *            find the branch multiplier with input tree with branch length
	 *            provided
	 * @param constFile
	 *            file name topological constraints
	 * @param quiet
	 *            runs PhyML in quiet mode
	 *
	 * @return String to call PhyML with all flags
	 */

	private static String formatCommand(boolean nesiSubmit, String seqFileName, String dataType,
			String sequential, String nDataSets, String parsimony,
			String bootstrap, String subModel, String aaRateFileName,
			String freq, String tstvRatio, String pInvar, String nSubCat,
			String gamma, String useMedian, String freeRates, String codePos,
			String search, String treeFile, String params, String randStart,
			String nRandStart, String randSeed, String printSitelnl,
			String printTrace, String runID, String noMemoryCheck,
			String noJcolAlias, String contrainedLens, String constFile,
			String quiet) {

		String command = "";
		if(!nesiSubmit){
			command += phymlVersion();
			if (command.equals("")) {
				return "";
			}
			if (!seqFileName.equals(""))
				command += "-i " + seqFileName + " ";
		}

		// if the string param is not "", add flag
		// this of course assumes that the things passed to it are ready to use
		// in
		// PhyML
		if (!dataType.equals(""))
			command += "-d " + dataType + " ";
		if (!sequential.equals(""))
			command += "-q ";
		if (!nDataSets.equals(""))
			command += "-n " + nDataSets + " ";
		if (!parsimony.equals(""))
			command += "-p ";
		if (!bootstrap.equals(""))
			command += "-b " + bootstrap + " ";
		if (!subModel.equals(""))
			command += "-m " + subModel + " ";
		if (!aaRateFileName.equals(""))
			command += "--aa_rate_file " + aaRateFileName + " ";
		if (!freq.equals(""))
			command += "-f " + freq + " ";
		if (!tstvRatio.equals(""))
			command += "-t " + tstvRatio + " ";
		if (!pInvar.equals(""))
			command += "-v " + pInvar + " ";
		if (!nSubCat.equals(""))
			command += "-c " + nSubCat + " ";
		if (!gamma.equals(""))
			command += "-a " + gamma + " ";
		if (!useMedian.equals("") && !useMedian.equals("mean"))
			command += "--use_median ";
		if (!freeRates.equals(""))
			command += "--free_rates ";
		if (!codePos.equals(""))
			command += "--codpos " + codePos + " ";
		if (!search.equals(""))
			command += "-s " + search + " ";
		if (!treeFile.equals(""))
			command += "-u " + treeFile + " ";
		if (!params.equals(""))
			command += "-o " + params + " ";
		if (!randStart.equals(""))
			command += "--rand_start ";
		if (!nRandStart.equals(""))
			command += "--n_rand_starts " + nRandStart + " ";
		if (!randSeed.equals(""))
			command += "--r_seed " + randSeed + " ";
		if (!printSitelnl.equals(""))
			command += "--print_site_lnl ";
		if (!printTrace.equals(""))
			command += "--print_trace ";
		if (!runID.equals(""))
			command += "--run_id " + runID + " ";
		if (!noMemoryCheck.equals(""))
			command += "--no_memory_check ";
		if (!noJcolAlias.equals(""))
			command += "--no_jcolalias ";
		if (!contrainedLens.equals(""))
			command += "--contrained_lens" + contrainedLens + " ";
		if (!constFile.equals(""))
			command += "--constraint_file " + constFile + " ";
		if (!quiet.equals(""))
			command += "--quiet ";

		return command;
	} // end formatCommand

	/**
	 * Method to obtain the correct path of the PhyML version to use depending
	 * on the OS and OS architecture of user computer. OS arch in practice is
	 * actually the bit version of the JVM.
	 *
	 * @return String of the absolute path of PhyML
	 */
	// Method to identify the OS and bit-version of OS.
	// returns the appropriate phyML version to execute
	private static String phymlVersion() {
		String phymlVer = "";

		// Identify the current OS to select correct version of PhyML to use
		String os = System.getProperty("os.name");
		System.out.println("OS: " + os);

		File binDir = new File(System.getProperty("user.dir"), "bin");

		myLogger.debug("Checking bin dir: "+binDir.getAbsolutePath());

		if ( ! binDir.exists() ) {
			binDir = new File(System.getProperty("user.dir")).getParentFile();
			binDir = new File(binDir, "phyml");
			myLogger.debug("Didn't work, now checking bin dir: "+binDir.getAbsolutePath());
			if ( ! binDir.exists() ) {
				throw new RuntimeException("Can't find directory for phyml binaries");
			}
		}

		myLogger.debug("Using bin dir: "+binDir.getAbsolutePath());

		// Windows OS
		if (os.matches("(?i).*Windows.*")) {
			phymlVer = binDir.getAbsolutePath()
					+ "\\PhyML-aBayes_win32.exe ";

			// command = phymlPath + " ";
		}
		// Mac OS
		if (os.matches("(?i).*Mac OS.*")) {
			// stick code in here to run OSX version of phyml
			phymlVer = binDir.getAbsolutePath()
					+ "/PhyML-aBayes_macOS_i386 ";
		}
		// Some Linux distro
		if (os.matches("(i?).*Linux.*")) {
			// see if 32 or 64 bit version
			String osArch = System.getProperty("os.arch");
			if (osArch.endsWith("64")) {
				System.out.println("Running 64-bit JVM, using 64-bit PhyML");
				phymlVer = binDir.getAbsolutePath()
						+ "/PhyML-Bayes_linux64 ";
				//phymlVer = "phyml ";
			} else {
				System.out.println("Running 32-bit JVM, using 32-bit PhyML");
				phymlVer = binDir
						+ "/PhyML-Bayes_linux32 ";
			}
		}
		if (phymlVer.equals("")) {
			System.out.println(os + " is not a supported operating system");
		}

		return phymlVer;
	} // end phymlVersion

} // end class
