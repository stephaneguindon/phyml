package phyml;

import java.awt.Color;

import javax.swing.JPanel;

/**
 * Implements all components to select the molecule type. This includes the
 * JPanels FileUpload, DataType, SeqFile and NumDataSets.
 * 
 * @author Christoph Knapp
 * 
 */
public class InputDataPanel extends JPanel {
	/**
	 * default id
	 */
	private static final long serialVersionUID = 1L;
	private FileUpload fU;
	private DataTypePhyml dT;
	private SeqFile sF;
	private NumDataSets nDS;

	/**
	 * Constructor initialises all objects and sets there sizes and positions.
	 */
	public InputDataPanel() {
		CustomGridLayout layout = new CustomGridLayout();
		setLayout(layout);
		layout.setDimensions(1, 0.15);
		add(new Separator("Input Data", Color.LIGHT_GRAY, null,false));
		layout.setDimensions(1, 0.15);
		fU = new FileUpload("");
		add(fU);
		layout.setDimensions(1, 0.02);
		add(new Separator(null, null, false));
		layout.setDimensions(1, 0.15);
		dT = new DataTypePhyml();
		add(dT);
		layout.setDimensions(1, 0.02);
		add(new Separator(null, null, false));
		layout.setDimensions(1, 0.15);
		sF = new SeqFile();
		add(sF);
		layout.setDimensions(1, 0.02);
		add(new Separator(null, null, false));
		layout.setDimensions(1, 0.30);
		nDS = new NumDataSets();
		add(nDS);
		layout.setDimensions(1, 0.04);
		add(new Separator(true, false));
	}

	/**
	 * Retrieves from the DataType object which molecule type is currently
	 * selected.
	 * 
	 * @return true if DNA, false if AA
	 */
	public Boolean isDNA() {
		return dT.isDNA();
	}

	/**
	 * Retrieves the type of molecule from the Datatype object.
	 * 
	 * @return String : "DNA" or "AA".
	 */
	public String getDataType() {
		return dT.getDataType();
	}

	/**
	 * Retrieves the file format from the SeqFile object.
	 * 
	 * @return "Interleaved" or "Sequential" depending on what is selected in
	 *         the gui.
	 */
	public String getFileFormat() {
		return sF.getFileFormat();
	}

	/**
	 * Retrieves the number of data sets in the input file as specified by the
	 * user from a NumDataSets object.
	 * 
	 * @return String representation of number of data sets.
	 */
	public String getNumDataSets() {
		return nDS.getNumDataSets();
	}

	/**
	 * Retrieves the run id from a NumDataSets object.
	 * 
	 * @return String : run id as typed in by user.
	 */
	public String getRunID() {
		return nDS.getRunID();
	}

	/**
	 * Retrives the file path to the input data set from a FileUpload object.
	 * 
	 * @return String : file path to input file.
	 */
	public String getInputPath() {
		return fU.getInputPath();
	}

	/**
	 * Forwards the request of changing the moleculetype to a DataType object.
	 * 
	 * @param b
	 *            boolean : true if DNA false if AA.
	 */
	public void SetDna(boolean b) {
		dT.SetDna(b);
		fU.setDna(b);
	}

	/**
	 * Forwards a request to a DataType object and states whether the components
	 * for selecting the molecule type are enabled or not.
	 * 
	 * @param b
	 *            boolean : true if components are enabled, false otherwise.
	 */
	public void UnLockDna(boolean b) {
		dT.UnLockDna(b);
	}

	/**
	 * Forwards the request of setting the file format to interleaved or
	 * sequential.
	 * 
	 * @param b
	 *            boolean : true if interleaved, false if sequential
	 */
	public void SetInterleaved(boolean b) {
		sF.SetInterleaved(b);
	}

	/**
	 * Forwards a request to a SeqFile object whether the interleaved/sequential
	 * tickboxes ar enabled or not.
	 * 
	 * @param b
	 *            boolean : true if enabled, false if disabled.
	 */
	public void UnLockInterleaved(boolean b) {
		sF.UnLockInterleaved(b);
	}

	/**
	 * Forwards a request to a NumDataSets object to reset the number of data
	 * sets back to 1.
	 */
	public void SetNumDS() {
		nDS.SetNumDS();
	}

	/**
	 * Forwards a request to a NumDataSets to enable or disable the components
	 * to change the number of data sets.
	 * 
	 * @param b
	 *            boolean : true if enabled, false if disabled
	 */
	public void UnLockNumDS(boolean b) {
		nDS.UnLockNumDS(b);
	}

	/**
	 * Requests from a FileUpload object whether the input file is user selected
	 * or the example data file.
	 * 
	 * @return boolean : true if user selected file, false if example file.
	 */
	public boolean getIsFile() {
		return fU.getIsFile();
	}

	/**
	 * Forwards a request to change the path in the JTextfield which specifies
	 * the input file.
	 * 
	 * @param inputPath
	 * 
	 *            String : file path to the input file.
	 */
	public void setInputPath(String inputPath) {
		fU.setInputPath(inputPath);
	}
}
