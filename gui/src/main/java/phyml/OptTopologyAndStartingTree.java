package phyml;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.DefaultComboBoxModel;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;

/**
 * Class for implementing all "Optimise Tree Topology" and "Starting tree"
 * components as well as setting their size and location.
 * 
 * @author Christoph Knapp
 * 
 */

public class OptTopologyAndStartingTree extends JPanel implements
		ActionListener {
	/**
	 * default id
	 */
	private static final long serialVersionUID = 1L;
//	private JComboBox optTopBox;
	private JComboBox staTreBox;
	private String[] staTreArr1;
	private String[] staTreArr2;
	private String userDefStartTreePath;
	private JRadioButton choice1;
	private JRadioButton choice2;

	/**
	 * Constructor to instatiate all components and setting their size and
	 * location.
	 */
	public OptTopologyAndStartingTree() {
		userDefStartTreePath = "";
		choice1 = new JRadioButton("Yes");
		choice2 = new JRadioButton("No");
		choice1.setSelected(true);
		choice1.addActionListener(this);
		choice2.addActionListener(this);
//		optTopBox = new JComboBox(new String[] { "yes", "no" });
//		optTopBox.addActionListener(this);
		staTreArr1 = new String[] { "BioNJ", "parsimony", "user tree" };
		staTreArr2 = new String[] { "BioNJ", "user tree" };
		staTreBox = new JComboBox(staTreArr1);
		staTreBox.addActionListener(this);
		
		CustomGridLayout layout = new CustomGridLayout();
		setLayout(layout);
		layout.setDimensions(1, 0.1);
		add(new JPanel());
		layout.setDimensions(0.01, 0.4);
		add(new JPanel());
		layout.setDimensions(0.33, 0.4);
		add(new JLabel("Optimise Tree Topology"));
		layout.setDimensions(0.27, 0.4);
		add(new JPanel());
		layout.setDimensions(0.38, 0.4);
		JPanel p1 = new JPanel();
		CustomGridLayout lo1 = new CustomGridLayout();
		p1.setLayout(lo1);
		add(p1);
		lo1.setDimensions(0.2, 1);
		p1.add(new JPanel());
		lo1.setDimensions(0.4, 1);
		p1.add(choice1);
		lo1.setDimensions(0.4, 1);
		p1.add(choice2);
		layout.setDimensions(0.01, 0.4);
		add(new JPanel());
		layout.setDimensions(1, 0.1);
		add(new JPanel());
		layout.setDimensions(0.01, 0.4);
		add(new JPanel());
		layout.setDimensions(0.33, 0.4);
		add(new JLabel("Starting Tree"));
		layout.setDimensions(0.27, 0.4);
		add(staTreBox);
		layout.setDimensions(0.39, 0.4);
		add(new JPanel());
//		
//		layout.setDimensions(0.44, 0.9);
//		JPanel p2 = new JPanel();
//		add(p2);
//		layout.setDimensions(0.01, 0.9);
//		add(new JPanel());
//		CustomGridLayout lO1 = new CustomGridLayout();
//		p1.setLayout(lO1);
//		lO1.setDimensions(0.62, 1);
//		p1.add(new JLabel("Optimise Tree Topology"));
//		lO1.setDimensions(0.1, 1);
//		p1.add(new JLabel());
//		lO1.setDimensions(0.28, 1);
//		p1.add(optTopBox);
//		CustomGridLayout lO2 = new CustomGridLayout();
//		p2.setLayout(lO2);
//		lO2.setDimensions(0.61, 1);
//		p2.add(new JLabel("Starting Tree"));
//		lO2.setDimensions(0.1, 1);
//		p2.add(new JLabel());
//		lO2.setDimensions(0.23, 1);
//		p2.add(staTreBox);
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		if(e.getSource() == choice1){
			int index = staTreBox.getSelectedIndex();
			if(choice1.isSelected()){
				choice2.setSelected(false);
				staTreBox.setModel(new DefaultComboBoxModel(staTreArr1));
				if (index == 1) {
					staTreBox.setSelectedIndex(2);
				}
				PhymlPanel.tS.setOptBraLenToYes(true);
			}else{
				choice2.setSelected(true);
				staTreBox.setModel(new DefaultComboBoxModel(staTreArr2));
				if (index == 1) {
					staTreBox.setSelectedIndex(0);
				} else if (index == 2) {
					staTreBox.setSelectedIndex(1);
				}
				PhymlPanel.tS.setOptBraLenToYes(false);
			}
		}else if(e.getSource() == choice2){
			int index = staTreBox.getSelectedIndex();
			if(choice2.isSelected()){
				choice1.setSelected(false);
				staTreBox.setModel(new DefaultComboBoxModel(staTreArr2));
				if (index == 1) {
					staTreBox.setSelectedIndex(0);
				} else if (index == 2) {
					staTreBox.setSelectedIndex(1);
				}
				PhymlPanel.tS.setOptBraLenToYes(false);
			}else{
				choice1.setSelected(true);
				staTreBox.setModel(new DefaultComboBoxModel(staTreArr1));
				if (index == 1) {
					staTreBox.setSelectedIndex(2);
				}
				PhymlPanel.tS.setOptBraLenToYes(true);
			}
		} else if (e.getSource() == staTreBox) {
			if (staTreBox.getSelectedItem().toString().equals("user tree")) {
				JFileChooser fc = new JFileChooser();
				fc.setCurrentDirectory(new File(System.getProperty("user.dir")));
				int returnVal = fc.showOpenDialog(this);
				if (returnVal == JFileChooser.APPROVE_OPTION) {
					userDefStartTreePath = fc.getSelectedFile()
							.getAbsolutePath();
				}
				if (userDefStartTreePath.equals("")) {
					staTreBox.setSelectedIndex(0);
				}
			} else {
				userDefStartTreePath = "";
			}
		}
	}

	/**
	 * What starting tree to use.
	 * 
	 * @return 
	 * String : "" if BioNJ or a user tree is selected, "parsimony"
	 * parsimony is selected.
	 */
	public String getStartTree() {
		if (staTreBox.getSelectedIndex() == 0
				|| staTreBox.getSelectedIndex() == 2) {
			return "";
		}
		return staTreBox.getSelectedItem().toString();
	}

	/**
	 * Retrieves the path to the user tree file.
	 * 
	 * @return 
	 * String : full path to the user defined starting tree file.
	 */
	public String getTreeFile() {
		if (staTreBox.getSelectedIndex() > 1) {
			return userDefStartTreePath;
		}
		return "";
	}

	/**
	 * Whether the tree topology should be optimised or not.
	 * 
	 * @return 
	 * boolean : true if the tree toplogy is supposed to be optimised,
	 * otherwise false.
	 */
	public boolean isOptimiseTreeTopology() {
		if (choice1.isSelected()) {
			return true;
		}
		return false;
	}
}
