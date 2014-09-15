package phyml;

import java.awt.Color;

import javax.swing.JPanel;
/**
 * Instantiats all components for the settings of tree
 * searching and to specify there size and location.
 * 
 * @author Christoph Knapp
 *
 */
public class TreeSearchingPanel extends JPanel {
	/**
	 * default id
	 */
	private static final long serialVersionUID = 1L;
	public OptTopologyAndStartingTree oTAST;
	private OptBraLenAndTreTop oBLATT;
	private AddRanStaTreAndNumRanStaTre aRSTANRST;

	/**
	 * Constructor for instantiating all components for the settings of tree
	 * searching and to specify there size and location.
	 */
	public TreeSearchingPanel() {
		CustomGridLayout layout = new CustomGridLayout();
		setLayout(layout);
		layout.setDimensions(1, 0.14);
		add(new Separator("Tree Searching", Color.lightGray, null,false));
		layout.setDimensions(1, 0.28);
		oTAST = new OptTopologyAndStartingTree();
		add(oTAST);
		layout.setDimensions(1, 0.28);
		oBLATT = new OptBraLenAndTreTop(true);
		add(oBLATT);
		layout.setDimensions(1, 0.28);
		aRSTANRST = new AddRanStaTreAndNumRanStaTre(true);
		add(aRSTANRST);
	}

	/**
	 * Retrieves the type of starting tree from an OptTopologyAndStartingTree
	 * object.
	 * 
	 * @return String : Empty String ("") if "BioNJ" or "user tree" is selected
	 *         and "parsimony" if parsimony is selected.
	 */
	public String getStartTree() {
		return oTAST.getStartTree();
	}

	/**
	 * Retrieves the selected search method for Tree Topology search.
	 * 
	 * @return String : "BEST" if "NNI and SPR" is selected and "SPR" if "SPR"
	 *         is selected, ""(empty String) otherwise.
	 */
	public String getSearch() {
		return oBLATT.getSearch();
	}

	/**
	 * Retrieves a file path to an user defined starting tree.
	 * 
	 * @return String : either the path to a file or "";
	 */
	public String getTreeFile() {
		return oTAST.getTreeFile();
	}

	/**
	 * Retrieves from an OptTopologyAndStartingTree object whether the tree
	 * topology should be optimised or not.
	 * 
	 * @return boolean : true if tree topoly should be optimised, otherwise
	 *         false.
	 */
	public boolean isOptimiseTreeTopology() {
		return oTAST.isOptimiseTreeTopology();
	}

	/**
	 * Retrieves from an OptBraLenAndTreTop object whether the branch length of
	 * a tree is optimised or not.
	 * 
	 * @return boolean : true if branch length should be optimised, otherwise
	 *         false.
	 */
	public boolean isOptimiseBranchLength() {
		return oBLATT.isOptimiseBranchLength();
	}

	/**
	 * Retrieves from an AddRanStaTreAndNumRanStaTre object whether a random
	 * starting tree is used or not.
	 * 
	 * @return String : Empty String for no random starting tree and yes if a
	 *         random starting tree is used.
	 */
	public String getRandStart() {
		return aRSTANRST.getRandStart();
	}

	/**
	 * Retrieves the number of random starting trees from an
	 * AddRanStaTreAndNumRanStaTre object.
	 * 
	 * @return String : Number of random starting tree as String.
	 */
	public String getNumRandStart() {
		return aRSTANRST.getNumRandStart();
	}

	/**
	 * Sets the ComboBox for "Opt Branch Length" to yes or no.
	 * 
	 * @param b
	 *            boolean : true if yes, false if no
	 */
	public void setOptBraLenToYes(boolean b) {
		oBLATT.setIsYes(b);
	}

	/**
	 * Sets the Combobox for "Tree Topology Search" to "NNI"
	 * 
	 * @param b
	 *            boolean : true if "NNI", otherwise false
	 */
	public void setIsNNI(boolean b) {
		aRSTANRST.setIsNNI(b);
	}
}
