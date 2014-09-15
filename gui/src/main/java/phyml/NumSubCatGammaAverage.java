package phyml;

import javax.swing.*;
import javax.swing.text.AttributeSet;
import javax.swing.text.BadLocationException;
import javax.swing.text.Document;
import javax.swing.text.PlainDocument;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
/**
 * Implements all components to specify the "# Subcategories",
 * "Gamma", "Alpha", and the "Average".
 *
 * @author Christoph Knapp
 *
 */
public class NumSubCatGammaAverage extends JPanel implements ActionListener,
		FocusListener {
	/**
	 * default id
	 */
	private static final long serialVersionUID = 1L;
	public CustomTextField subCat;
//	private JComboBox gammaBox;
	private CustomTextField alpha;
//	private JComboBox average;
	private JLabel lab2;
	private JLabel lab1;
	private JLabel lab3;
	private JLabel lab4;
	private JRadioButton gammaSahpeEstimated;
	private JRadioButton gammaShapeFixed;
	private JRadioButton averageMeanRadioButton;
	private JRadioButton averageMedianRadioButton;

	/**
	 * Constructor method implements all components to specify the "# Subcategories", "Gamma",
	 * "Alpha", and the "Average", it instantiates all components, adds listeners sets size
	 * and position of components and adds the components.
	 */
	public NumSubCatGammaAverage() {
		subCat = new CustomTextField("", false);
        subCat.setText("4");
//		gammaBox = new JComboBox(new String[] { "estimated", "fixed" });
//		gammaBox.addActionListener(this);
		gammaSahpeEstimated = new JRadioButton("Estimated");
		gammaSahpeEstimated.setSelected(true);
		gammaShapeFixed = new JRadioButton("Fixed");
		gammaSahpeEstimated.addActionListener(this);
		gammaShapeFixed.addActionListener(this);
		alpha = new CustomTextField("", true);
        alpha.setText("1.0");
		alpha.setEnabled(false);
		alpha.addFocusListener(this);
//		average = new JComboBox(new String[] { "mean", "median" });
		averageMeanRadioButton = new JRadioButton("Mean");
		averageMeanRadioButton.setSelected(true);
		averageMedianRadioButton = new JRadioButton("Median");
		averageMeanRadioButton.addActionListener(this);
		averageMedianRadioButton.addActionListener(this);
		lab1 = new JLabel("Number of Rate Classes");
		lab2 = new JLabel("Gamma Shape Parameter");
		lab3 = new JLabel("Value");
		lab4 = new JLabel("Average");

		CustomGridLayout layout = new CustomGridLayout();
		setLayout(layout);
		layout.setDimensions(1, 0.1);
		add(new JPanel());
		layout.setDimensions(0.01, 0.24);
		add(new JPanel());
		layout.setDimensions(0.33, 0.24);
		add(lab1);
		layout.setDimensions(0.27, 0.24);
		add(subCat);
		layout.setDimensions(0.39, 0.24);
		add(new JPanel());
		layout.setDimensions(1, 0.1);
		add(new JPanel());
		layout.setDimensions(0.01, 0.24);
		add(new JPanel());
		layout.setDimensions(0.33, 0.24);
		add(lab2);
		layout.setDimensions(0.27, 0.24);
		add(alpha);
		layout.setDimensions(0.38, 0.24);
		JPanel p1 = new JPanel();
		CustomGridLayout lo1 = new CustomGridLayout();
		p1.setLayout(lo1);
		add(p1);
		lo1.setDimensions(0.2, 1);
		p1.add(new JPanel());
		lo1.setDimensions(0.4, 1);
		p1.add(gammaSahpeEstimated);
		lo1.setDimensions(0.4, 1);
		p1.add(gammaShapeFixed);
		layout.setDimensions(0.01, 0.24);
		add(new JPanel());
		layout.setDimensions(1, 0.1);
		add(new JPanel());
		layout.setDimensions(0.01, 0.24);
		add(new JPanel());
		layout.setDimensions(0.33, 0.24);
		add(lab4);
		layout.setDimensions(0.27, 0.24);
		add(new JPanel());
		layout.setDimensions(0.38, 0.24);
		JPanel p2 = new JPanel();
		CustomGridLayout lo2 = new CustomGridLayout();
		p2.setLayout(lo2);
		add(p2);
		lo2.setDimensions(0.2, 1);
		p2.add(new JPanel());
		lo2.setDimensions(0.4, 1);
		p2.add(averageMeanRadioButton);
		lo2.setDimensions(0.4, 1);
		p2.add(averageMedianRadioButton);
		layout.setDimensions(0.01, 0.24);
		add(new JPanel());

	}

	/**
	 * Sets all Components either to visible or unvisible.
	 *
	 * @param b
	 *            boolean : true if visible, false otherwise.
	 */
	public void setCompVisible(boolean b) {
		subCat.setVisible(b);
		gammaSahpeEstimated.setVisible(b);
		gammaShapeFixed.setVisible(b);
		averageMeanRadioButton.setVisible(b);
		averageMedianRadioButton.setVisible(b);
		alpha.setVisible(b);
		lab1.setVisible(b);
		lab2.setVisible(b);
		lab3.setVisible(b);
		lab4.setVisible(b);
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		if(e.getSource() == gammaSahpeEstimated){
			if(gammaSahpeEstimated.isSelected()){
				gammaShapeFixed.setSelected(false);
				alpha.setEnabled(false);
			}else{
				gammaShapeFixed.setSelected(true);
				alpha.setEnabled(true);
				alpha.requestFocus();
			}
		}
		if(e.getSource() == gammaShapeFixed){
			if(gammaShapeFixed.isSelected()){
				gammaSahpeEstimated.setSelected(false);
				alpha.setEnabled(true);
				alpha.requestFocus();
			}else{
				gammaSahpeEstimated.setSelected(true);
				alpha.setEnabled(false);
			}
		}
		if(e.getSource() == averageMeanRadioButton){
			if(averageMeanRadioButton.isSelected()){
				averageMedianRadioButton.setSelected(false);
			}else{
				averageMedianRadioButton.setSelected(true);
			}
		}
		if(e.getSource() == averageMedianRadioButton){
			if(averageMedianRadioButton.isSelected()){
				averageMeanRadioButton.setSelected(false);
			}else{
				averageMeanRadioButton.setSelected(true);
			}
		}
	}

    public void disableGammaShape(boolean disabled) {
            gammaSahpeEstimated.setEnabled(!disabled);
            gammaShapeFixed.setEnabled(!disabled);
    }

    public void setGammaShapeFixed(boolean fixed) {
        gammaSahpeEstimated.setSelected(!fixed);
        gammaShapeFixed.setSelected(fixed);
    }

    class CustomTextField extends JTextField {
		/**
		 * default id
		 */
		private static final long serialVersionUID = 1L;
		private boolean type;

		/**
		 * Constructor method.
		 *
		 * @param string
		 *            String : double value to set as default in text field.
		 * @param type
		 *            boolean : whether it is an integer or double. True if
		 *            Double value, false if not.
		 */
		public CustomTextField(String string, boolean type) {
			super(string);
			this.type = type;
			if (type) {
				if (!isDouble(string)) {
					this.setText("");
				}
			} else {
				if (!isInteger(string)) {
					this.setText("");
				}
			}
		}

		/**
		 * Tests whether the input String is convertable into an Integer.
		 *
		 * @param str
		 *            String : input String.
		 * @return
		 * boolean : true if convertible, false otherwise.
		 */
		private boolean isInteger(String str) {
			try {
				Integer.parseInt(str);
				return true;
			} catch (NumberFormatException e) {
				return false;
			}
		}

		@Override
		protected Document createDefaultModel() {
			return new doubleOnlyDocument();
		}

		/**
		 * Tests whether the input String is convertable into an Double.
		 *
		 * @param str
		 *            String : input String.
		 * @return
		 * boolean : true if convertible, false otherwise.
		 */
		private boolean isDouble(String str) {
			try {
				Double.parseDouble(str);
				return true;
			} catch (NumberFormatException e) {
				if (str.equals(".") && !getGammaString().contains(".")) {
					return true;
				}
				return false;
			}
		}

		/**
		 *
		 * @author Christoph Knapp
		 *
		 */
		class doubleOnlyDocument extends PlainDocument {
			/**
			 * default id
			 */
			private static final long serialVersionUID = 1L;

			@Override
			public void insertString(int offs, String str, AttributeSet a)
					throws BadLocationException {
				if (type) {
					if (str == null || !isDouble(str)) {
						return;
					}
					super.insertString(offs, str, a);
				} else {
					if (str == null || !isInteger(str)) {
						return;
					}
					super.insertString(offs, str, a);
				}
			}
		}
	}

	/**
	 * Retrieves the String of the Gamma value from the Costume TextField.
	 *
	 * @return
	 * String : Gama value as a String.
	 */
	public String getGammaString() {
		return alpha.getText();
	}

	/**
	 * Retrieves the number of subCategories as a String.
	 *
	 * @return
	 * String : "4" by default otherwise what the user typed in.
	 */
	public String getNumSubCat() {
		if (subCat.getText().equals("")) {
			return 4 + "";
		}
		return subCat.getText();
	}

	/**
	 * Retrieves the Alpha value set by the user.
	 *
	 * @return
	 * String : returns "e" for estimated and the Alpha value as
	 * specified by the user otherwise.
	 */
	public String getAlpha() {
		if (gammaSahpeEstimated.isSelected()) {
			return "e";
		} else if (gammaShapeFixed.isSelected()) {
			return alpha.getText();
		}
		return "e";
	}

	/**
	 * Retrieves what average type is selected.
	 *
	 * @return
	 * String : "mean" or "median".
	 */
	public String getUseMedian() {
		String ret = "";
		if(averageMeanRadioButton.isSelected()){
			ret="mean";
		}else if(averageMedianRadioButton.isSelected()){
			ret="median";
		}
		return ret;
	}

	@Override
	public void focusGained(FocusEvent e) {
	}

	@Override
	public void focusLost(FocusEvent e) {
		if (alpha.getText().equals("")) {
			gammaSahpeEstimated.setSelected(true);
			gammaShapeFixed.setSelected(false);
			alpha.setEnabled(false);
		}
	}
}
