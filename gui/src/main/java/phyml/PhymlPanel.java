package phyml;

import grisu.control.ServiceInterface;
import grisu.frontend.control.login.LoginException;
import grisu.frontend.control.login.LoginManager;
import grith.gridsession.GridClient;
import grith.jgrith.cred.AbstractCred;
import grith.jgrith.cred.Cred;
import nz.org.nesi.phyml.swing.GridPanel;
import phyml.view.CredCreationDialog;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

/**
 *
 * This Class implements the main JPanel . It specifies the layout of the
 * different components and controls the interaction between the different
 * components.
 *
 * @author Christoph Knapp
 * @version 15-June-2012
 *
 */

public class PhymlPanel extends JPanel implements ActionListener {

    private CredCreationDialog ccd = null;

    public void showLoginDialog() {
        getCredCreationDialog().setVisible(true);
    }

    private CredCreationDialog getCredCreationDialog() {
        if ( ccd == null ) {
            ccd = new CredCreationDialog(gridClient);
            ccd.pack();
        }
        return ccd;
    }

	/**
	 * default id
	 */
	private static final long serialVersionUID = 1L;
	public static InputDataPanel iDP;
	public static SubstitutionModelPhyml sM;
	public static TreeSearchingPanel tS;
	private BranchSupport bS;
	private static JButton submit;
	private StandardOutPanel standardOut;
	private GridPanel gridPanel;
//	private static TreePanel treePanel;
	private JTabbedPane tabbedPane;
	private JCheckBox nesiSubmit;
    private final GridClient gridClient;

	/**
	 * This Class is the main panel for the graphical user interface. I
	 * specifies the sub panels the phyml gui is subdivided in.
	 */
	public PhymlPanel(GridClient gc) {

        gridClient = gc;
        CustomGridLayout mainLayout = new CustomGridLayout();
		setLayout(mainLayout);
		mainLayout.setDimensions(1, 1);
		tabbedPane = new JTabbedPane();
		JPanel mainPan = new JPanel();
		submit = new JButton("Submit");
		submit.addActionListener(this);
		CustomGridLayout layout = new CustomGridLayout();
		mainPan.setLayout(layout);
		iDP = new InputDataPanel();
		layout.setDimensions(1, 0.2);
		mainPan.add(iDP);
		if (iDP.isDNA()) {
			sM = new SubstitutionModelPhyml("DNA");
		} else {
			sM = new SubstitutionModelPhyml("AA");
		}
        // need to do that here because apparently it should be the default
        // any other position results in npe
        sM.setNumSubCatGammaAverageCompVisible(false);
		layout.setDimensions(1, 0.37);
		mainPan.add(sM);
		tS = new TreeSearchingPanel();
		layout.setDimensions(1, 0.23);
		mainPan.add(tS);
		bS = new BranchSupport();
		layout.setDimensions(1, 0.13);
		mainPan.add(bS);
		layout.setDimensions(1, 0.07);
		JPanel p1 = new JPanel();
		mainPan.add(p1);
		CustomGridLayout lO1 = new CustomGridLayout();
		p1.setLayout(lO1);
		lO1.setDimensions(1, 0.2);
		p1.add(new JPanel());
		lO1.setDimensions(0.35, 0.6);
		p1.add(new Separator(false));
		lO1.setDimensions(0.3, 0.6);
		p1.add(submit);
		lO1.setDimensions(0.1, 0.6);
		p1.add(new JPanel());
		lO1.setDimensions(0.25, 0.6);
		nesiSubmit = new JCheckBox("Submit to NeSI");
		p1.add(nesiSubmit);
		//p1.add(new Separator(false));
		lO1.setDimensions(1, 0.2);
		p1.add(new JPanel());
		layout.setDimensions(1, 0.01);
		mainPan.add(new Separator(false));
		gridPanel = new GridPanel(gridClient);
		tabbedPane.addTab("Run PhyML", mainPan);
		standardOut = new StandardOutPanel();
		tabbedPane.addTab("Standard output", standardOut);
		tabbedPane.addTab("NeSI", gridPanel);

//        tabbedPane.addChangeListener(new ChangeListener() {
//            public void stateChanged(ChangeEvent e) {
////                System.out.println("Tab: " + tabbedPane.getSelectedIndex());
//                if ( tabbedPane.getSelectedIndex() == 2 ) {
//                    if ( gridPanel.getServiceInterface() != null ) {
//                        RunningJobManager.getDefault(gridPanel.getServiceInterface()).updateJobnameList("PhyML", true);
//                    }
//                }
//            }
//        });
		//settings =  new SettingsPanel();
		//tabbedPane.addTab("Settings", new JPanel());
//		treePanel = new TreePanel();
//		tabbedPane.addTab("Trees", treePanel);
		add(tabbedPane);
	}

	@Override
	public void actionPerformed(ActionEvent a) {

        if ( nesiSubmit.isSelected() ) {

            if (gridPanel.getServiceInterface() == null) {

                Cred c = AbstractCred.getExistingCredential();

                if (c == null || !c.isValid()) {
                    getCredCreationDialog().setVisible(true);

                }

                ServiceInterface si = null;
                try {
                    si = LoginManager.login("bestgrid", gridClient.getCredential(), false);
                } catch (LoginException e) {
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }

                gridPanel.setServiceInterface(si);
            }

        }

		if (!iDP.getInputPath().equals("") && (new File(iDP.getInputPath())).exists()) {
//			if(!nesiSubmit.isSelected()){
				tabbedPane.setSelectedIndex(1);
//			}
			callPHYML.callPhyML(gridPanel, nesiSubmit.isSelected(),
				iDP.getInputPath(),
				/**
				 * -i (or --input) seq file name seq file name is the
				 * name of the nucleotide or amino-acid sequence file in
				 * PHYLIP format.
				 */
				iDP.getDataType(),
				/**
				 * -d (or --datatype) data type data type is nt for
				 * nucleotide (default) and aa for amino-acid sequences.
				 */
				iDP.getFileFormat(),
				/**
				 * -q (or --sequential) Changes interleaved format
				 * (default) to sequential format.
				 */
				iDP.getNumDataSets(),
				/**
				 * -n (or --multiple) nb data sets nb data sets is an
				 * integer giving the number of data sets to analyse.
				 */
				tS.getStartTree(),// TODO check this method
				/**
				 * -p (or --pars)Use a minimum parsimony starting tree.
				 * This option is taken into account when the ’-u’
				 * option is absent and when tree topology modifications
				 * are to be done.
				 */
				bS.getNumBootstrap(),// TODO check this method
				/**
				 * -b (or --bootstrap) int – int > 0: int is the number
				 * of bootstrap replicates. – int = 0: neither
				 * approximate likelihood ratio test nor bootstrap
				 * values are computed. – int = -1: approximate
				 * likelihood ratio test returning aLRT statistics. –
				 * int = -2: approximate likelihood ratio test returning
				 * Chi2-based para-metric branch supports. – int = -4:
				 * SH-like branch supports alone.
				 */
				sM.getSubModel(),
				/**
				 * -m (or --model) model name model name : substitution
				 * model name. – Nucleotide-based models: HKY85
				 * (default) | JC69 | K80 | F81 | F84 | TN93 | GTR |
				 * custom The custom option can be used to define a new
				 * substitution model. A string of six digits identifies
				 * the model. For instance, 000000 corresponds to F81
				 * (or JC69 provided the distribution of nucleotide
				 * frequencies is uni-form). 012345 corresponds to GTR.
				 * This option can be used for encoding any model that
				 * is a nested within GTR. See Section 6.1.2. NOTE: the
				 * substitution parameters of the custom model will be
				 * optimised so as to maximise the likelihood. It is
				 * possible to specify and fix (i.e., avoid
				 * optimisation) the values of the substitution rates
				 * only through the PHYLIP-like interface. – Amino-acid
				 * based models: LG (default) WAG | JTT | MtREV |
				 * Dayhoff | DCMut | RtREV | CpREV | VT | Blosum62 |
				 * MtMam | MtArt | HIVw | HIVb | custom The custom
				 * option is useful when one wants to use an amino-acid
				 * substitution model that is not available by default
				 * in PhyML. The symmetric part of the rate matrix, as
				 * well as the equilibrium amino-acid frequencies, are
				 * given in a file which name is asked for by the
				 * program. The format of this file is described in
				 * section 7.4.
				 */
				sM.getAARateFileName(),
				/**
				 * --aa rate file file name This option is compulsory
				 * when analysing amino-acid sequences under a ‘custom’
				 * model. file name should provide a rate matrix and
				 * equilibrium amino acid in PAML format (see Section ).
				 */
				sM.getFrequencies(),
				/**
				 * -f e, m, or “fA,fC,fG,fT” Nucleotide or amino-acid
				 * frequencies. – e : the character frequencies are
				 * determined as follows : - Nucleotide sequences:
				 * (Empirical) the equilibrium base frequencies are
				 * estimated by counting the occurence of the different
				 * bases in the alignment. - Amino-acid sequences:
				 * (Empirical) the equilibrium amino-acid frequencies
				 * are estimated by counting the occurence of the
				 * different amino-acids in the alignment. – m : the
				 * character frequencies are determined as follows : -
				 * Nucleotide sequences: (ML) the equilibrium base
				 * frequencies are estimated using maximum likelihood. -
				 * Amino-acid sequences: (Model) the equilibrium
				 * amino-acid frequencies are estimated using the
				 * frequencies defined by the substitution model. –
				 * “fA,fC,fG,fT” : only valid for nucleotide-based
				 * models. fA, fC, fG and fT are floating numbers that
				 * correspond to the frequencies of A, C, G and T
				 * respectively.
				 */
				sM.getTsTvRatio(),
				/**
				 * -t (or --ts/tv) ts/tv ratio ts/tv ratio:
				 * transition/transversion ratio. DNA sequences only.
				 * Can be a fixed positive value (e.g., 4.0) or type e
				 * to get the maximum likelihood estimate.
				 */
				sM.getPropInvarSites(),
				/**
				 * -v (or --pinv) prop invar prop invar: proportion of
				 * invariable sites. Can be a fixed value in the [0,1]
				 * range or type e to get the maximum likelihood
				 * estimate.
				 */
				sM.getNumSubCat(),
				/**
				 * -c (or --nclasses) nb subst cat nb subst cat: number
				 * of relative substitution rate categories. Default: nb
				 * subst cat=4. Must be a positive integer.
				 */
				sM.getAlpha(),
				/**
				 * -a (or --alpha) gamma gamma: value of the gamma shape
				 * parameter. Can be a fixed positive value or e to get
				 * the maximum likelihood estimate. The value of this
				 * parameter is estimated in the maximum likelihood
				 * framework by default.
				 */
				sM.getUseMedian(),
				/**
				 * --use median The middle of each substitution rate
				 * class in the discrete gamma distribution is taken as
				 * the median. The mean is used by default.
				 */
				"",// TODO check with Stephane
				/**
				 * --free rates As an alternative to the discrete gamma
				 * model, it is possible to estimate the (relative) rate
				 * in each class of the (mixture) model and the
				 * corresponding frequencies. This model has more
				 * parameters than the discrete gamma one but usually
				 * provides a significantly better fit to the data.
				 */
				"",// TODO check with Stephane
				/**
				 * --codpos 1,2 or 3 When analysing an alignment of
				 * coding sequences, use this option to consider only
				 * the first, second or third coding position for the
				 * estimation.
				 */
				tS.getSearch(),
				/**
				 * -s (or --search) move Tree topology search operation
				 * option. Can be either NNI (default, fast) or SPR (a
				 * bit slower than NNI) or BEST (best of NNI and SPR
				 * search).
				 */
				tS.getTreeFile(),// TODO check this method
				/**
				 * -u (or --inputtree) user tree file user tree file:
				 * starting tree filename. The tree must be in Newick
				 * format.
				 */
				getParams(),
				/**
				 * -o params This option focuses on specific parameter
				 * optimisation. – params=tlr: tree topology (t), branch
				 * length (l) and substitution rate parameters (r) are
				 * optimised. – params=tl: tree topology and branch
				 * lengths are optimised. – params=lr: branch lengths
				 * and substitution rate parameters are optimised. –
				 * params=l: branch lengths are optimised. – params=r:
				 * substitution rate parameters are optimised. –
				 * params=n: no parameter is optimised.
				 */
				tS.getRandStart(),
				/**
				 * --rand start This option sets the initial tree to
				 * random. It is only valid if SPR searches are to be
				 * performed.
				 */
				tS.getNumRandStart(),
				/**
				 * --n rand starts num num is the number of initial
				 * random trees to be used. It is only valid if SPR
				 * searches are to be performed.
				 */
				"",// TODO randSeed (might generate a random number
					// here)
				/**
				 * --r seed num num is the seed used to initiate the
				 * random number generator. Must be an integer.
				 */
				"yes",// TODO printSitelnl (for now set to yes)
				/**
				 * --print site lnl Print the likelihood for each site
				 * in file * phyml lk.txt.
				 */
				"yes",// TODO printTrace (for now set to yes)
				/**
				 * --print trace Print each phylogeny explored during
				 * the tree search process in file * phyml trace.txt.
				 */
				iDP.getRunID(),
				/**
				 * --run id ID string Append the string ID string at the
				 * end of each PhyML output file. This optionmay be
				 * useful when running simulations involving PhyML. It
				 * can also be used to ‘tag’ multiple analysis of the
				 * same data set with various program settings.
				 */
				"yes",// TODO noMemoryCheck (for now set to yes later on
						// we could add a warning and let the user
						// decide)
				/**
				 * --no memory check By default, when processing a large
				 * data set, PhyML will pause and ask the user to
				 * confirm that she/he wants to continue with the
				 * execution of the analysis despite the large amount of
				 * memory required. The --no memory check skips this
				 * question. It is especially useful when running PhyML
				 * in batch mode.
				 */
				"",// TODO noJcolAlias (for now we want it faster)
				/**
				 * -no jcolalias By default, PhyML preprocesses each
				 * alignment by putting together (or aliasing) the
				 * columns that are identical. Use this option to skip
				 * this step but be aware that the analysis might then
				 * take more time to complete.
				 */
				"",// TODO contrainedLens for now empty need to check
					// with Stephane
				/**
				 * --contrained lens When an input tree with branch
				 * lengths is provided, this option will find the branch
				 * multiplier that maximises the likelihood (i.e., the
				 * relative branch lengths remain constant)
				 */
				"",// TODO constFile not provided yet
				/**
				 * --constraint file file name file name lists the
				 * topological constraints under which the tree topology
				 * search is conducted. This option should be used in
				 * conjunction with -u file name. See Section 7.5 in the
				 * phyml manual for more information.
				 */
				""// TODO quiet
			);
			submit.setEnabled(false);
		}
	}

	/**
	 * Checks which parameters are selected in the gui and forms the input
	 * String for the "-o params" parameter
	 *
	 * @return – params=tlr: tree topology (t), branch length (l) and
	 *         substitution rate parameters (r) are optimised.<br>
	 *         – params=tl: tree topology and branch lengths are optimised.<br>
	 *         – params=lr: branch lengths and substitution rate parameters are
	 *         optimised.<br>
	 *         – params=l: branch lengths are optimised.<br>
	 *         – params=r: substitution rate parameters are optimised.<br>
	 *         – params=n: no parameter is optimised.<br>
	 */
	private String getParams() {
		if (tS.isOptimiseTreeTopology() && tS.isOptimiseBranchLength()
				&& sM.isOptimisedRateParameter()) {
			return "tlr";
		} else if (tS.isOptimiseTreeTopology() && tS.isOptimiseBranchLength()) {
			return "tl";
		} else if (tS.isOptimiseBranchLength() && sM.isOptimisedRateParameter()) {
			return "lr";
		} else if (tS.isOptimiseBranchLength()) {
			return "l";
		} else if (sM.isOptimisedRateParameter()) {
			return "r";
		} else if (!tS.isOptimiseTreeTopology() && !tS.isOptimiseBranchLength()
				&& !sM.isOptimisedRateParameter()) {
			return "n";
		}
		return "";
	}

	/**
	 * Loads and refreshes the tree panel this method is called when a datafile
	 * is specified and when an analysis finishes.
	 *
	 */
	public static void loadTrees() {
		File f = new File(iDP.getInputPath());
//		treePanel.setFiles(f.getParent());
	}

	/**
	 * Disables or enables the submit button. If an analysis is running it is
	 * disabled otherwise not
	 *
	 */
	public static void SetSubmit(boolean b) {
		submit.setEnabled(true);
	}

	/**
	 * Sets the molecule type in the InputDataPanel object
	 *
	 * @param b
	 *            : true if DNA false if AA
	 */
	public static void SetDna(boolean b) {
		iDP.SetDna(b);
	}

	/**
	 * Sets the tickboxes for "Data Type" to disabled if the example file is
	 * selected and enabled if not
	 *
	 * @param b
	 *            : true if a file is uploaded, false if a the "example file"
	 *            tickbox is selected.
	 */
	public static void UnLockDna(boolean b) {
		iDP.UnLockDna(b);
	}

	/**
	 * Sets the tickboxes for "Sequence File"
	 *
	 * @param b
	 *            true if File is in interleaved ansd false if the File is in
	 *            sequential
	 */
	public static void SetInterleaved(boolean b) {
		iDP.SetInterleaved(b);
	}

	/**
	 * Enables or disables the "Sequence File" tickboxes
	 *
	 * @param b
	 *            : true if enabled, false if disabled
	 */
	public static void UnLockInterleaved(boolean b) {
		iDP.UnLockInterleaved(b);
	}

	/**
	 * Sets the number of datasets back to 1, this method is called when the
	 * example data file is used.
	 */
	public static void SetNumDS() {
		iDP.SetNumDS();
	}

	/**
	 * Enables or disables the JSpinner for the number of datasets.
	 *
	 * @param b
	 *            true if enabled, false if disabled
	 */
	public static void UnLockNumDS(boolean b) {
		iDP.UnLockNumDS(b);
	}

	/**
	 * Sets the JTextField for the input path. this is used when the user
	 * chooses a file through the JFileChooser or the "Example File" tickbox is
	 * selected.
	 *
	 * @param inputPath
	 */
	public static void setInputFile(String inputPath) {
		iDP.setInputPath(inputPath);
	}
	/**
	 * Passes on the right moleculetype to all the subclasses i.e substitutionModel
	 *
	 * @param molType
	 * String : either "DNA" or "AA"
	 */
	public static void setMoleculeType(String molType) {
		sM.setMoleculeType(molType);
		if(molType.equals("DNA")){
			SetDna(true);
		}else{
			SetDna(false);
		}
	}
	public static void setRatioBoxOff(boolean off){
		sM.setRatioBoxOff(off);
	}

	public static void setRatioBoxDefault() {
		sM.setRatioBoxDefault();
	}
}
