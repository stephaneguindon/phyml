package phyml;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;

/**
 * Extends JPanel and implements all components for specifying the format of the
 * file.
 * 
 * @author Christoph Knapp
 * 
 */

public class SeqFile extends JPanel implements ActionListener {
	/**
	 * default id
	 */
	private static final long serialVersionUID = 1L;
	private JLabel lab1;
	private JRadioButton choice1;
	private JRadioButton choice2;

	/**
	 * Constructor method implements all components for specifying the format of
	 * the file.
	 */
	public SeqFile() {
		lab1 = new JLabel("Sequence File Format");
		choice1 = new JRadioButton("Interleaved");
		choice1.setToolTipText("If selected, the input file is in Interleaved format.");
		choice1.setSelected(true);
		choice2 = new JRadioButton("Sequential");
		choice2.setToolTipText("If selected, the input file is in Sequential format.");
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
			} else {
				choice2.setSelected(true);
			}
		}
		if (e.getSource() == choice2) {
			if (choice2.isSelected()) {
				choice1.setSelected(false);
			} else {
				choice1.setSelected(true);
			}
		}
	}

	/**
	 * Whether "Interleaved" or "Sequential" is selected.
	 * 
	 * @return boolean : true if "Interleaved", false if "Sequential"
	 */
	public boolean isInterleaved() {
		return this.choice1.isSelected();
	}

	/**
	 * Retrieves the file format as a String.
	 * 
	 * @return String : "Interleaved" or "Sequential";
	 */
	public String getFileFormat() {
		if (choice1.isSelected()) {
			return "Interleaved";
		} else {
			return "Sequential";
		}
	}

	/**
	 * Changes the file format.
	 * 
	 * @param b
	 *            boolean : true if "Interleaved", false for "Sequential"
	 */
	public void SetInterleaved(boolean b) {
		if (b) {
			choice1.setSelected(true);
			choice2.setSelected(false);
		} else {
			choice1.setSelected(false);
			choice2.setSelected(true);
		}
	}

	/**
	 * Enables or disables the components for selecting a file format.
	 * 
	 * @param b
	 *            boolean : true if enabled, false otherwise
	 */
	public void UnLockInterleaved(boolean b) {
		choice1.setEnabled(b);
		choice2.setEnabled(b);
	}
}
