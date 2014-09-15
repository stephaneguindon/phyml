package phyml;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.JTextField;
import javax.swing.SpinnerModel;
import javax.swing.SpinnerNumberModel;

/**
 * Class to instantiate the JSpinner for specifying the number of data sets and
 * the Run ID.
 * 
 * @author Christoph Knapp
 * 
 */

public class NumDataSets extends JPanel {
	/**
	 * default id
	 */
	private static final long serialVersionUID = 1L;
	private JLabel lab1;
	private JLabel lab2;
	private JTextField iDTxt;
	private JSpinner spinner;

	/**
	 * Constructor method instantiate the JSpinner for specifying the number of 
	 * data sets and the Run ID.
	 */
	public NumDataSets() {
		lab1 = new JLabel();
		lab1.setText("Number Of Data Sets");
		SpinnerModel model = new SpinnerNumberModel(1,// initial value
				1,// min
				null,// max
				1);// step
		spinner = new JSpinner(model);
		lab2 = new JLabel("Run ID");
		iDTxt = new JTextField();
		
		CustomGridLayout layout = new CustomGridLayout();
		setLayout(layout);
		layout.setDimensions(1, 0.1);
		add(new JPanel());
		layout.setDimensions(0.01, 0.4);
		add(new JPanel());
		layout.setDimensions(0.33, 0.4);
		add(lab1);
		layout.setDimensions(0.27, 0.4);
		add(spinner);
		layout.setDimensions(0.38, 0.4);
		add(new JPanel());
		layout.setDimensions(0.01, 0.4);
		add(new JPanel());
		layout.setDimensions(1, 0.1);
		add(new JPanel());
		layout.setDimensions(0.01, 0.4);
		add(new JPanel());
		layout.setDimensions(0.33, 0.4);
		add(lab2);
		layout.setDimensions(0.27, 0.4);
		add(iDTxt);
		layout.setDimensions(0.38, 0.4);
		add(new JPanel());
	}

	/**
	 * Retrieves the number of data sets as specified by the user.
	 * 
	 * @return 
	 * String : The String representation of the selected Integer in the JSpinner.
	 */
	public String getNumDataSets() {
		return ""
				+ ((SpinnerNumberModel) spinner.getModel()).getNumber()
						.intValue();
	}

	/**
	 * Retrieves the text typed in to the JTextField for "Run ID"
	 * 
	 * @return 
	 * String : Text typed into "Run ID" JTextField.
	 */
	public String getRunID() {
		return iDTxt.getText();
	}

	/**
	 * Resets the number in the JSpinner to 1.
	 */
	public void SetNumDS() {
		spinner.setValue(1);
	}

	/**
	 * Enables or disables the JSpinner for specifying the number of data sets.
	 * 
	 * @param b
	 * boolean : If true the JSpinner is enabled otherwise its disabled.
	 */
	public void UnLockNumDS(boolean b) {
		spinner.setEnabled(b);
	}
}
