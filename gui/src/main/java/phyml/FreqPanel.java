package phyml;

import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.text.AttributeSet;
import javax.swing.text.BadLocationException;
import javax.swing.text.Document;
import javax.swing.text.PlainDocument;

/**
 * Implements the JPanel for specifying the frequencies of DNA nucleotides.
 * 
 * @author Christoph Knapp
 * 
 */

/**
 * @author christoph
 *
 */
/**
 * @author christoph
 *
 */
public class FreqPanel extends JPanel {
	/**
	 * default id
	 */
	private static final long serialVersionUID = 1L;
	private CustomTextField t1;
	private CustomTextField t2;
	private CustomTextField t3;
	private CustomTextField t4;
	private boolean hasAllValues;
	private JLabel l1;

	/**
	 * Instantiates the Components displayed in the graphical user interface for
	 * setting the frequencies of DNA nucleotides.
	 */
	public FreqPanel() {
		hasAllValues = false;
		t1 = new CustomTextField();
		t2 = new CustomTextField();
		t3 = new CustomTextField();
		t4 = new CustomTextField();
		l1 = new JLabel("Freq A - Freq C - Freq G - Freq T");
		CustomGridLayout layout = new CustomGridLayout();
		setLayout(layout);
		
		layout.setDimensions(0.01, 1);
		add(new JPanel());
		layout.setDimensions(0.33, 1);
		add(l1);
		JPanel p1 = new JPanel();
		CustomGridLayout lo1 = new CustomGridLayout();
		p1.setLayout(lo1);
		layout.setDimensions(0.27, 1);
		add(p1);
		lo1.setDimensions(0.45, 1);
		p1.add(t1);
		lo1.setDimensions(0.1, 1);
		p1.add(new JPanel());
		lo1.setDimensions(0.45, 1);
		p1.add(t2);
		
		
		JPanel p2 = new JPanel();
		CustomGridLayout lo2 = new CustomGridLayout();
		p2.setLayout(lo2);
		layout.setDimensions(0.38, 0.9);
		add(p2);
		lo2.setDimensions(0.2, 1);
		p2.add(new JPanel());
		lo2.setDimensions(0.35, 1);
		p2.add(t3);
		lo2.setDimensions(0.08, 1);
		p2.add(new JPanel());
		lo2.setDimensions(0.33, 1);
		p2.add(t4);
	}

	/**
	 * Instantiates the Components displayed in the graphical user interface for
	 * setting the frequencies of DNA nucleotides.
	 * 
	 * @param d1
	 *            : Sets the default value displayed in the JTextField.
	 * @param d2
	 *            : Sets the default value displayed in the JTextField.
	 * @param d3
	 *            : Sets the default value displayed in the JTextField.
	 * @param d4
	 *            : Sets the default value displayed in the JTextField.
	 */
	public FreqPanel(double d1, double d2, double d3, double d4) {
		t1 = new CustomTextField();
		t2 = new CustomTextField();
		t3 = new CustomTextField();
		t4 = new CustomTextField();
		if ((d1 + d2 + d3 + d4) == 1) {
			t1.setText("" + d1);
			t2.setText("" + d1);
			t3.setText("" + d1);
			t4.setText("" + d1);
			hasAllValues = true;
		}
		l1 = new JLabel("Freq A - Freq C - Freq G - Freq T");
		CustomGridLayout layout = new CustomGridLayout();
		setLayout(layout);
		
		layout.setDimensions(0.01, 1);
		add(new JPanel());
		layout.setDimensions(0.33, 1);
		add(l1);
		
		JPanel p1 = new JPanel();
		CustomGridLayout lo1 = new CustomGridLayout();
		p1.setLayout(lo1);
		layout.setDimensions(0.27, 1);
		add(p1);
		lo1.setDimensions(0.45, 1);
		p1.add(t1);
		lo1.setDimensions(0.1, 1);
		add(new JPanel());
		layout.setDimensions(0.45, 1);
		p1.add(t2);
		JPanel p2 = new JPanel();
		CustomGridLayout lo2 = new CustomGridLayout();
		p2.setLayout(lo2);
		layout.setDimensions(0.38, 1);
		add(p1);
		lo2.setDimensions(0.2, 1);
		p2.add(new JPanel());
		lo2.setDimensions(0.35, 1);
		p2.add(t3);
		layout.setDimensions(0.1, 1);
		add(new JPanel());
		layout.setDimensions(0.35, 1);
		p2.add(t4);
	}

	/**
	 * Tests whether the input variable in the JTextFields are valid.
	 * 
	 * @return true if the sum of all typed in values is <=1 and false if not.
	 */
	public boolean isValidSet() {
		double d1 = convertToDouble(t1.getText());
		double d2 = convertToDouble(t2.getText());
		double d3 = convertToDouble(t3.getText());
		double d4 = convertToDouble(t4.getText());
		if (hasAllValues) {
			if ((d1 + d2 + d3 + d4) == 1) {
				return true;
			}
		} else {
			if ((d1 + d2 + d3 + d4) <= 1) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Class, implementing a custom JTextField which only makes it possible to
	 * type in values smaller than 1 and greater than 0.
	 * 
	 * @author Christoph Knapp
	 * 
	 */
	private class CustomTextField extends JTextField implements FocusListener {
		/**
		 * default id
		 */
		private static final long serialVersionUID = 1L;

		/**
		 * Constructor method only adds FocusListener.
		 */
		public CustomTextField() {
			this.addFocusListener(this);
		}
		@Override
		protected Document createDefaultModel() {
			return new doubleOnlyDocument();
		}

		/**
		 * Extends the PlainDocument class, makes it possible to only type in
		 * Double values smaller or equal 1 but not less than 0.
		 */
		class doubleOnlyDocument extends PlainDocument {
			/**
			 * default id
			 */
			private static final long serialVersionUID = 1L;

			@Override
			public void insertString(int offs, String str, AttributeSet a)
					throws BadLocationException {
				if (str == null || !isNumber(str) || !isRightInput(offs, str)) {
					return;
				}
				if ((str.contains(".") && !this.getText(0, this.getLength())
						.contains(".")) || !str.contains(".")) {
					super.insertString(offs, str, a);
				}
			}

			/**
			 * Tests if the input so far is valid
			 * 
			 * @param str
			 *            : What is typed into the text field while its typed
			 *            in.
			 * @return true if input is valid false if not.
			 */
			private boolean isNumber(String str) {
				try {
					if(str.length()>1){
						Double.parseDouble(str);
						return true;
					}
					Integer.parseInt(str);
					return true;
				} catch (NumberFormatException e) {
					if (str.equals(".")) {
						return true;
					}
					return false;
				}
			}
		}

		@Override
		public void focusGained(FocusEvent e) {
		}

		/**
		 * Tests how the String starts and if the it is possible that it will be
		 * a valid double value.
		 * 
		 * @param offs
		 *            : position of typed segment.
		 * @param str
		 *            : segment which was just typed in.
		 * @return true if input is valid ,false if not
		 */
		public boolean isRightInput(int offs, String str) {
			String txt = this.getText();
			if (offs == 0
					&& (str.equals("0") || str.equals(".") || str.equals("1"))) {
				return true;
			} else if (offs == 1 && txt.charAt(0) == '.') {
				return true;
			} else if (offs == 1 && txt.charAt(0) == '0' && str.equals(".")) {
				return true;
			} else if (offs > 1) {
				return true;
			} else if (offs > 0 && txt.charAt(0) == '1') {
				return false;
			}else if(str.length()>1){
				if(grep(str, '.')>1){
					return false;
				}
				for(char c : str.toCharArray()){
					try{
						if(c!='.'){
							Integer.parseInt(c+"");
						}
					}catch(NumberFormatException e){
						return false;
					}
				}
				return true;
			}
			return false;
		}

		private int grep(String str, char match) {
			int count = 0;
			for(char c : str.toCharArray()){
				if(c==match){
					count++;
				}
			}
			return count;
		}
		@Override
		public void focusLost(FocusEvent e) {
			if (!isDouble(this.getText())) {
				this.setText("");
				hasAllValues = false;
			}
			if (checkHasAllValues()) {
				hasAllValues = true;
				if (!isValidSet()) {
					hasAllValues = false;
					this.setText("");
				}
			}
		}
	}

	/**
	 * Tests if all values typed into the text fields are valid.
	 * 
	 * @return true if they are valid false otherwise.
	 */
	public boolean checkHasAllValues() {
		return isDouble(t1.getText()) && isDouble(t2.getText())
				&& isDouble(t3.getText()) && isDouble(t4.getText());
	}

	/**
	 * Converts the a string to a double value or 0 if it is not possible.
	 * 
	 * @param text
	 *            : String to be converted into Double.
	 * @return The converted Double value or 0 if not possible to convert.
	 */
	private double convertToDouble(String text) {
		try {
			return Double.parseDouble(text);
		} catch (NumberFormatException e) {
			return 0;
		}
	}

	/**
	 * Tests whether a String is convertible to Double.
	 * 
	 * @param text
	 *            : String to be tested.
	 * @return True if possible, false if not.
	 */
	public boolean isDouble(String text) {
		try {
			Double.parseDouble(text);
			return true;
		} catch (NumberFormatException e) {
			return false;
		}
	}

	/**
	 * Enables or disables all Frequency text fields.
	 * 
	 * @param b
	 *            : If true components are enabled, if not disabled.
	 */
	public void setCompEnabled(boolean b) {
		t1.setEnabled(b);
		t2.setEnabled(b);
		t3.setEnabled(b);
		t4.setEnabled(b);
	}

	/**
	 * Sets all Text Fields including Labels visible or not.
	 * 
	 * @param b
	 *            : If true all components are visible, otherwise not.
	 */
	public void setCompVisible(boolean b) {
		l1.setVisible(b);
		t1.setVisible(b);
		t2.setVisible(b);
		t3.setVisible(b);
		t4.setVisible(b);
	}

	/**
	 * Retrieves the input from the Text Fields and converts it into the right
	 * format.
	 * 
	 * @return : String array containing 4 elements, one for every Frequency. If
	 *         not possible this method returns null.
	 */
	public String[] getFreq() {
		if (hasAllValues) {
			return new String[] { Double.parseDouble(t1.getText()) + "",
					Double.parseDouble(t4.getText()) + "",
					Double.parseDouble(t4.getText()) + "",
					Double.parseDouble(t4.getText()) + "" };
		} else if (isValidSet()) {
			double[] values = new double[] { convertToDouble(t1.getText()),
					convertToDouble(t2.getText()),
					convertToDouble(t3.getText()),
					convertToDouble(t4.getText()) };
			boolean[] whichZeroes = new boolean[] { false, false, false, false };
			if ((values[0] + values[1] + values[2] + values[3]) < 1) {
				int count = 0;
				if (values[0] == 0) {
					whichZeroes[0] = true;
					count++;
				}
				if (values[1] == 0) {
					whichZeroes[1] = true;
					count++;
				}
				if (values[2] == 0) {
					whichZeroes[2] = true;
					count++;
				}
				if (values[3] == 0) {
					whichZeroes[3] = true;
					count++;
				}
				double value = (1 - values[0] + values[1] + values[2] + values[3])
						/ count;
				for (int i = 0; i < whichZeroes.length; i++) {
					if (whichZeroes[i]) {
						values[i] = value;
					}
				}
				return new String[] { "" + values[0], "" + values[1],
						"" + values[2], "" + values[3] };
			} else {
				return new String[] { "" + values[0], "" + values[1],
						"" + values[2], "" + values[3] };
			}
		}
		return null;
	}
	/**
	 * Sets the frequencies of the frequency text fields to the parameters passed in.
	 * @param freqA
	 * 				: Frequency for the A nucleotide
	 * @param freqC
	 * 				: Frequency for the C nucleotide
	 * @param freqG
	 * 				: Frequency for the G nucleotide
	 * @param freqT
	 * 				: Frequency for the T nucleotide
	 */
	public void setFrequencies(String freqA, String freqC, String freqG,String freqT) {
		t1.setText(freqA);
		t2.setText(freqC);
		t3.setText(freqG);
		t4.setText(freqT);
	}
}
