package phyml;

import grisu.jcommons.utils.PackageFileHelper;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.io.File;

import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTextField;

/**
 * Implements the components to upload a file or to select an example file.
 * 
 * @author Christoph Knapp
 */
public class FileUpload extends JPanel implements ActionListener, FocusListener {
	/**
	 * default id
	 */
	private static final long serialVersionUID = 1L;
	private String exampleFilePath;
	private JTextField pathField;
	private JButton browse;
	private boolean hasChoice;
	private JRadioButton choice1;
	private JRadioButton choice2;

	/**
	 * Constructor for instantiating all components and set there position and
	 * size.
	 * 
	 * @param exampleFilePath
	 *            String : file path to an input file.
	 */
	public FileUpload(String exampleFilePath) {
		this.exampleFilePath = exampleFilePath;
		choice1 = new JRadioButton("File");
		choice1.setToolTipText("If selected, the choosen input file is used.");
		choice1.setSelected(true);
		choice2 = new JRadioButton("Example File");
		choice2.setToolTipText("If selected, the example input file is used.");
		choice2.setSelected(false);
		pathField = new JTextField(exampleFilePath);
		browse = new JButton("Browse");
		hasChoice = false;
		JLabel lab1;
		lab1 = new JLabel("Sequences (Phylip Format)");
		
		CustomGridLayout layout = new CustomGridLayout();
		setLayout(layout);
		layout.setDimensions(1, 0.1);
		add(new JPanel());
		layout.setDimensions(0.01, 0.9);
		add(new JPanel());
		layout.setDimensions(0.33, 0.9);
		add(lab1);
		JPanel p1 = new JPanel();
		layout.setDimensions(0.27, 0.9);
		add(p1);
		CustomGridLayout lO1 = new CustomGridLayout();
		p1.setLayout(lO1);
		lO1.setDimensions(0.49, 1);
		p1.add(pathField);
		lO1.setDimensions(0.03, 1);
		p1.add(new JPanel());
		lO1.setDimensions(0.45, 1);
		p1.add(browse);
		lO1.setDimensions(0.03, 1);
		p1.add(new JPanel());
		layout.setDimensions(0.38, 0.9);
		JPanel p2 = new JPanel();
		add(p2);
		CustomGridLayout lO2 = new CustomGridLayout();
		p2.setLayout(lO2);
		lO2.setDimensions(0.2, 1);
		p2.add(new JPanel());
		lO2.setDimensions(0.4, 1);
		p2.add(choice1);
		lO2.setDimensions(0.4, 1);
		p2.add(choice2);
		layout.setDimensions(0.01, 0.9);
		add(new JPanel());
		layout.setDimensions(1, 0.05);
		add(new JPanel());
		browse.addActionListener(this);
		choice1.addActionListener(this);
		choice2.addActionListener(this);
		pathField.addFocusListener(this);
	}

	/**
	 * Returns the file path of the input data file.
	 * 
	 * @return String : path to input file.
	 */
	public String getInputPath() {
		if (!hasChoice) {
			return pathField.getText();
		} else {
			if (choice1.isSelected()) {
				return pathField.getText();
			} else {
				return exampleFilePath;
			}
		}
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		if (e.getSource() == choice1) {
			if (choice1.isSelected()) {
				pathField.setEnabled(true);
				browse.setEnabled(true);
				choice2.setSelected(false);
				PhymlPanel.UnLockDna(true);
				System.out.println("unlock");
				PhymlPanel.UnLockInterleaved(true);
				PhymlPanel.UnLockNumDS(true);
				PhymlPanel.setInputFile("");
			} else {
				pathField.setEnabled(false);
				browse.setEnabled(false);
				choice2.setSelected(true);
				PhymlPanel.setMoleculeType("DNA");
				PhymlPanel.UnLockDna(false);
				System.out.println("lock");
				PhymlPanel.SetInterleaved(true);
				PhymlPanel.UnLockInterleaved(false);
				PhymlPanel.SetNumDS();
				PhymlPanel.UnLockNumDS(false);
				
				String temp = PackageFileHelper.getPath("phyml_ex.txt");
				PhymlPanel.setInputFile(temp);
//				PhymlPanel.setInputFile(System.getProperty("user.dir")
//						+ "/exampleOutput/phyml_ex.txt");
				PhymlPanel.loadTrees();
			}
		} else if (e.getSource() == choice2) {
			if (choice2.isSelected()) {
				pathField.setEnabled(false);
				browse.setEnabled(false);
				choice1.setSelected(false);
				PhymlPanel.setMoleculeType("DNA");
				PhymlPanel.UnLockDna(false);
				System.out.println("lock");
				PhymlPanel.SetInterleaved(true);
				PhymlPanel.UnLockInterleaved(false);
				PhymlPanel.SetNumDS();
				PhymlPanel.UnLockNumDS(false);
				String temp = PackageFileHelper.getPath("phyml_ex.txt");
				PhymlPanel.setInputFile(temp);
//				PhymlPanel.setInputFile(System.getProperty("user.dir")
//						+ "/exampleOutput/phyml_ex.txt");
				PhymlPanel.loadTrees();
			} else {
				pathField.setEnabled(true);
				browse.setEnabled(true);
				choice1.setSelected(true);
				PhymlPanel.UnLockDna(true);
				System.out.println("unlock");
				PhymlPanel.UnLockInterleaved(true);
				PhymlPanel.UnLockNumDS(true);
				PhymlPanel.setInputFile("");
			}
		} else if (e.getSource() == browse) {
			JFileChooser fc = new JFileChooser();
			fc.setCurrentDirectory(new File(System.getProperty("user.dir")));
			int returnVal = fc.showOpenDialog(this);
			if (returnVal == JFileChooser.APPROVE_OPTION) {
				pathField.setText(fc.getSelectedFile().getAbsolutePath());
			}
			if (!pathField.getText().equals("")) {
				PhymlPanel.loadTrees();
			}
		}
	}

	@Override
	public void focusGained(FocusEvent e) {
	}

	@Override
	public void focusLost(FocusEvent e) {
		if (!pathField.getText().equals("")) {
			PhymlPanel.loadTrees();
		}
	}

	/**
	 * Returns whether an input data file is selected or the example file is
	 * used.
	 * 
	 * @return boolean : true if user defined file, false if example input file.
	 */
	public boolean getIsFile() {
		return !choice2.isSelected();
	}

	/**
	 * Changes the text inside the JTextField which stores the file path to the
	 * input data file.
	 * 
	 * @param inputPath
	 *            String : file path where the input data is located. Note, this
	 *            method does not check if the file is valid.
	 */
	public void setInputPath(String inputPath) {
		pathField.setText(inputPath);
	}

	public void setDna(boolean b) {
		choice1.setEnabled(b);
		choice2.setEnabled(b);
	}
}
