package phyml;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

/**
 * Implements all components necessary for
 * "Tree topology search and "Optimise branch length".
 *
 * @author Christoph Knapp
 */

public class OptBraLenAndTreTop extends JPanel implements ActionListener {
	/**
	 * default id
	 */
	private static final long serialVersionUID = 1L;
	private JComboBox TreTopBox;
	private JLabel lab1;
	private JLabel lab2;
	private JRadioButton choice1;
	private JRadioButton choice2;

	/**
	 * Implements all components necessary for
	 * "Tree topology search and "Optimise branch length".
	 *
	 * @param isYes
	 *            Whether the DropDown menu for selecting Optimise Branch length
	 *            is visible (false) or not (true).
	 */
	public OptBraLenAndTreTop(boolean isYes) {
		lab1 = new JLabel("Tree Topology Search");
		lab2 = new JLabel("Optimise Branch Length");
		TreTopBox = new JComboBox(new String[] { "NNI", "SPR", "NNI and SPR" });
		choice1 = new JRadioButton("Yes");
		choice2 = new JRadioButton("No");
		choice1.addActionListener(this);
		choice2.addActionListener(this);
		choice1.setSelected(true);
		TreTopBox.addActionListener(this);
		if (isYes) {
			choice1.setEnabled(false);
			choice2.setEnabled(false);
		} else {
			TreTopBox.setEnabled(false);
		}

		CustomGridLayout layout = new CustomGridLayout();
		setLayout(layout);
		layout.setDimensions(1, 0.1);
		add(new JPanel());
		layout.setDimensions(0.01, 0.4);
		add(new JPanel());
		layout.setDimensions(0.33, 0.4);
		add(lab1);
		layout.setDimensions(0.27, 0.4);
		add(TreTopBox);
		layout.setDimensions(0.39, 0.4);
		add(new JPanel());
		layout.setDimensions(1, 0.1);
		add(new JPanel());
		layout.setDimensions(0.01, 0.4);
		add(new JPanel());
		layout.setDimensions(0.33, 0.4);
		add(lab2);
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
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		if (TreTopBox.getSelectedItem().toString().equals("NNI")) {
			PhymlPanel.tS.setIsNNI(true);
		} else {
			PhymlPanel.tS.setIsNNI(false);
		}
		if(e.getSource() == choice1){
			if(choice1.isSelected()){
				choice2.setSelected(false);
			}else{
				choice2.setSelected(true);
			}
		}
		if(e.getSource() == choice2){
			if(choice2.isSelected()){
				choice1.setSelected(false);
			}else{
				choice1.setSelected(true);
			}
		}
	}

	/**
	 * Sets whether the Optimise Branch length components or the Tree topology
	 * components are visible.
	 *
	 * @param b
	 *            If true the optimise branch length components are not visible.
	 *            and the Tree topology components are visible.<br>
	 *            If false its the opposite from above.
	 */
	public void setIsYes(boolean b) {
		if (b) {
			choice1.setEnabled(false);
			choice2.setEnabled(false);
			TreTopBox.setEnabled(true);
			if (TreTopBox.getSelectedIndex() > 0) {
				PhymlPanel.tS.setIsNNI(false);
			} else {
				PhymlPanel.tS.setIsNNI(true);
			}
		} else {
			choice1.setEnabled(true);
			choice2.setEnabled(true);
			TreTopBox.setEnabled(false);
			PhymlPanel.tS.setIsNNI(true);
		}
	}

	/**
	 * Passes the user choice from the "Tree Topology Search" options.
	 *
	 * @return "BEST" if "NNI and SPR" is selected and "SPR" if "SPR" is
	 *         selected, ""(empty String) otherwise.
	 */
	public String getSearch() {
		if (PhymlPanel.tS.isOptimiseTreeTopology()) {
			if (TreTopBox.getSelectedItem().toString().equals("NNI and SPR")) {
				return "BEST";
			} else if (TreTopBox.getSelectedItem().toString().equals("SPR")) {
				return TreTopBox.getSelectedItem().toString();
			}
		}
		return "";
	}

	/**
	 * Retrieves whether the "Optimise Branch Length" parameter is set to yes or
	 * no.
	 *
	 * @return Returns true if "Optimise Branch Length" is set to yes otherwise
	 *         false.
	 */
	public boolean isOptimiseBranchLength() {
		if (choice1.isSelected()) {
			return true;
		}
		return false;
	}

}
