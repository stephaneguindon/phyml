package phyml;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;

/**
 * Class to implement all components for specifying the Molecule type "DNA" or
 * "AA".
 * 
 * @author Christoph Knapp
 * 
 */
public class DataTypePhyml extends JPanel implements ActionListener {
	/**
	 * default id
	 */
	private static final long serialVersionUID = 1L;
	private JLabel lab1;
	private JRadioButton choice1;
	private JRadioButton choice2;

	/**
	 * Constructor instantiates all components for specifying the Molecule type
	 * "DNA" or "AA".
	 */
	public DataTypePhyml() {
		lab1 = new JLabel("Data Type");
		choice1 = new JRadioButton("DNA");
		choice1.setToolTipText("If selected, the input type is specified as DNA.");
		choice1.setSelected(true);
		choice2 = new JRadioButton("Amino-acids");
		choice2.setToolTipText("If selected, the input type is specified as Amino-acids.");
		choice1.addActionListener(this);
		choice2.addActionListener(this);
		
		CustomGridLayout layout = new CustomGridLayout();
		setLayout(layout);
		layout.setDimensions(1, 0.1);
		add(new JPanel());
		layout.setDimensions(0.01, 0.9);
		add(new JPanel());
		layout.setDimensions(0.33, 0.9);
		add(lab1);
		layout.setDimensions(0.27, 0.9);
		add(new JPanel());
		layout.setDimensions(0.38, 0.9);
		JPanel p = new JPanel();
		add(p);
		CustomGridLayout lO1 = new CustomGridLayout();
		p.setLayout(lO1);
		lO1.setDimensions(0.2, 1);
		p.add(new JPanel());
		lO1.setDimensions(0.4, 1);
		p.add(choice1);
		lO1.setDimensions(0.4, 1);
		p.add(choice2);
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		if (e.getSource() == choice1) {
			if (choice1.isSelected()) {
				choice2.setSelected(false);
				PhymlPanel.setMoleculeType("DNA");
			} else {
				choice2.setSelected(true);
				PhymlPanel.setMoleculeType("AA");
			}
		}
		if (e.getSource() == choice2) {
			if (choice2.isSelected()) {
				choice1.setSelected(false);
				PhymlPanel.setMoleculeType("AA");
			} else {
				choice1.setSelected(true);
				PhymlPanel.setMoleculeType("DNA");
			}
		}
	}

	/**
	 * Tests whether the molecule type is set to DNA or not.
	 * 
	 * @return true if "DNA", false otherwise
	 */
	public boolean isDNA() {
		return this.choice1.isSelected();
	}

	/**
	 * retrieves the command line option for molecule type.
	 * 
	 * @return nt for DNA or aa for Amino Acids.
	 */
	public String getDataType() {
		if (choice1.isSelected()) {
			return "nt";
		} else {
			return "aa";
		}
	}

	/**
	 * Sets the components of the graphical user interface specifying the
	 * molecule type either to DNA or to AA (Amino Acid).
	 * 
	 * @param b
	 * boolean : if true the components are set to DNA otherwise AA.
	 */
	public void SetDna(boolean b) {
		if (b) {
			choice1.setSelected(true);
			choice2.setSelected(false);
		} else {
			choice1.setSelected(false);
			choice2.setSelected(true);
		}
	}

	/**
	 * Sets the components for choosing the molecule type to enabled or disabled
	 * 
	 * @param b
	 * boolean : If true components are enable otherwise disabled.
	 */
	public void UnLockDna(boolean b) {
		choice1.setEnabled(b);
		choice2.setEnabled(b);
	}
}
