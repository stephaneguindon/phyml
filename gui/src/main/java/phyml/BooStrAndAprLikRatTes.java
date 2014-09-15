package phyml;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;

import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTextField;
import javax.swing.text.AttributeSet;
import javax.swing.text.BadLocationException;
import javax.swing.text.Document;
import javax.swing.text.PlainDocument;

/**
 * Class for implementing all components to set the values for whether
 * "Non parametric Bootstrapping" is used, how many "Replicates", whether
 * additional statistics are printed to file and which
 * "Approximate Likelihood Ratio Test" is used.
 * 
 * @author Christoph Knapp
 * 
 */
public class BooStrAndAprLikRatTes extends JPanel implements ActionListener,
		FocusListener {
	/**
	 * default id
	 */
	private static final long serialVersionUID = 1L;
	private JComboBox likRatBox;
	private CustomTextField numRepField;
	private JRadioButton choice1;
	private JRadioButton choice2;
	private JRadioButton choice3;
	private JRadioButton choice4;

	/**
	 * Constructor method for implementing all components to set the values for
	 * whether "Non parametric Bootstrapping" is used, how many "Replicates",
	 * whether additional statistics are printed to file and which
	 * "Approximate Likelihood Ratio Test" is used.
	 */
	public BooStrAndAprLikRatTes() {
		choice1 = new JRadioButton("no");
		choice2 = new JRadioButton("yes");
		choice1.addActionListener(this);
		choice2.addActionListener(this);
		choice1.setSelected(true);
		likRatBox = new JComboBox(new String[] { "no", "SH-like supports",
				"aLRT statistics", "Chi2-based supports" });
		choice3 = new JRadioButton("no");
		choice4 = new JRadioButton("yes");
		choice3.addActionListener(this);
		choice4.addActionListener(this);
		choice3.setSelected(true);
		likRatBox.setSelectedIndex(1);
		numRepField = new CustomTextField("");
		numRepField.setEnabled(false);
		numRepField.addFocusListener(this);
		choice3.setEnabled(false);
		choice4.setEnabled(false);
		likRatBox.addActionListener(this);
		
		CustomGridLayout layout = new CustomGridLayout();
		setLayout(layout);
		layout.setDimensions(1, 0.1);
		add(new JPanel());
		layout.setDimensions(0.01, 0.24);
		add(new JPanel());
		layout.setDimensions(0.33, 0.24);
		add(new JLabel("Non Parametric Bootstrapping"));
		layout.setDimensions(0.27, 0.24);
		add(numRepField);
		layout.setDimensions(0.38, 0.24);
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
		layout.setDimensions(0.01, 0.24);
		add(new JPanel());
		layout.setDimensions(1, 0.1);
		add(new JPanel());
		layout.setDimensions(0.01, 0.24);
		add(new JPanel());
		layout.setDimensions(0.33, 0.24);
		add(new JLabel("Print Trees"));
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
		p2.add(choice3);
		lo2.setDimensions(0.4, 1);
		p2.add(choice4);
		layout.setDimensions(0.01, 0.24);
		add(new JPanel());
		layout.setDimensions(1, 0.1);
		add(new JPanel());
		layout.setDimensions(0.01, 0.24);
		add(new JPanel());
		layout.setDimensions(0.33, 0.24);
		add(new JLabel("Approximate Likelihood Ratio Test"));
		layout.setDimensions(0.27, 0.24);
		add(likRatBox);
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		if(e.getSource() == choice1){
			if(choice1.isSelected()){
				choice2.setSelected(false);
				numRepField.setEnabled(false);
				choice1.setEnabled(false);
				choice2.setEnabled(false);
			}else{
				choice2.setSelected(true);
				likRatBox.setSelectedIndex(0);
				numRepField.setEnabled(true);
				numRepField.requestFocus();
				choice1.setEnabled(true);
				choice2.setEnabled(true);
			}
		}
		if(e.getSource() == choice2){
			if(choice2.isSelected()){
				choice1.setSelected(false);
				likRatBox.setSelectedIndex(0);
				numRepField.setEnabled(true);
				numRepField.requestFocus();
				choice1.setEnabled(true);
				choice2.setEnabled(true);
			}else{
				choice1.setSelected(true);
				numRepField.setEnabled(false);
				choice1.setEnabled(false);
				choice2.setEnabled(false);
			}
		}
		if(e.getSource() == choice3){
			if(choice3.isSelected()){
				choice4.setSelected(false);
			}else{
				choice4.setSelected(true);
			}
		}
		if(e.getSource() == choice4){
			if(choice4.isSelected()){
				choice3.setSelected(false);
			}else{
				choice3.setSelected(true);
			}
		}
		if (e.getSource() == likRatBox) {
			if (likRatBox.getSelectedIndex() > 0) {
				numRepField.setEnabled(false);
				choice3.setEnabled(false);
				choice4.setEnabled(false);
				choice1.setSelected(true);
				choice2.setSelected(false);
			}
		}
	}

	/**
	 * Class to customize a JTextField that the user is only able to type in
	 * integer values.
	 * 
	 * @author Christoph Knapp
	 * 
	 */
	class CustomTextField extends JTextField {
		/**
		 * default id
		 */
		private static final long serialVersionUID = 1L;

		/**
		 * Class to customize a JTextField that the user is only able to type in
		 * integer values.
		 * 
		 * @param string
		 *            String : Default String as displayed when initialised.
		 */
		public CustomTextField(String string) {
			super(string);
			if (!isInteger(string)) {
				this.setText("");
			}
		}

		/**
		 * Tests whether the input value is convertible to an Integer value.
		 * 
		 * @param str
		 *            String : input String
		 * @return boolean : true if convertable, false otherwise.
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
		 * Class which extends PlainDocument for checking the input in the
		 * JTextField as the user types it in.
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
				if (str == null || !isInteger(str)) {
					return;
				}
				super.insertString(offs, str, a);
			}
		}
	}

	/**
	 * Retrieves the number of replicates if user selected to use Bootstrapping
	 * and typed in a number in the text field.
	 * 
	 * @return String : number of replicates.
	 */
	public String getNumBootstrap() {
		if (choice2.isSelected() && !numRepField.getText().equals("")) {
			return numRepField.getText();
		}
		return "";
	}

	@Override
	public void focusGained(FocusEvent e) {
	}

	@Override
	public void focusLost(FocusEvent e) {
		System.out.println("lost");
		if (choice2.isSelected()&&numRepField.getText().equals("")) {
			choice1.setSelected(true);
			choice2.setSelected(false);
			choice3.setSelected(true);
			choice4.setSelected(false);
			likRatBox.setSelectedIndex(1);
		}
	}
}
