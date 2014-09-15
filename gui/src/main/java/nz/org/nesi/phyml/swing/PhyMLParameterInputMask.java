package nz.org.nesi.phyml.swing;

import grisu.control.ServiceInterface;
import grisu.control.exceptions.RemoteFileSystemException;
import grisu.frontend.view.swing.jobcreation.widgets.Cpus;
import grisu.frontend.view.swing.jobcreation.widgets.SingleInputGridFile;
import grisu.frontend.view.swing.jobcreation.widgets.TextCombo;
import grisu.frontend.view.swing.jobcreation.widgets.Walltime;
import grisu.model.dto.GridFile;

import java.awt.Color;

import javax.swing.JPanel;
import javax.swing.border.LineBorder;
import javax.swing.border.TitledBorder;

import com.jgoodies.forms.layout.ColumnSpec;
import com.jgoodies.forms.layout.FormLayout;
import com.jgoodies.forms.layout.FormSpecs;
import com.jgoodies.forms.layout.RowSpec;

public class PhyMLParameterInputMask extends JPanel {
	
	private TextCombo textCombo;
	private SingleInputGridFile singleInputGridFile;
	private Cpus cpus;
	private Walltime walltime;
	private String command;
	private String file;

	public PhyMLParameterInputMask() {
		command=SwingClient.getCommand();
		file=SwingClient.getFile();
		System.out.println("PhyMLParameterInputMask "+command);
		System.out.println("PhyMLParameterInputMask "+file);
		setLayout(new FormLayout(new ColumnSpec[] {
				FormSpecs.RELATED_GAP_COLSPEC,
				ColumnSpec.decode("default:grow"),
				FormSpecs.RELATED_GAP_COLSPEC,
				FormSpecs.DEFAULT_COLSPEC,
				FormSpecs.RELATED_GAP_COLSPEC,},
			new RowSpec[] {
				FormSpecs.RELATED_GAP_ROWSPEC,
				FormSpecs.DEFAULT_ROWSPEC,
				FormSpecs.RELATED_GAP_ROWSPEC,
				FormSpecs.DEFAULT_ROWSPEC,
				FormSpecs.RELATED_GAP_ROWSPEC,
				FormSpecs.DEFAULT_ROWSPEC,}));
		add(getTextCombo(), "2, 2, 3, 1, fill, fill");
		add(getSingleInputGridFile(), "2, 4, 3, 1, fill, fill");
		add(getWalltime(), "2, 6, right, fill");
		add(getCpus(), "4, 6, fill, fill");
	}

	private TextCombo getTextCombo() {
		if (textCombo == null) {
			textCombo = new TextCombo();
			textCombo.setText(command);
			textCombo.setHistoryKey("PHYML_PARAMETERS");
			textCombo.setBorder(new TitledBorder(null, "phyML parameters", TitledBorder.LEADING, TitledBorder.TOP, null, null));
		}
		return textCombo;
	}
	private SingleInputGridFile getSingleInputGridFile() {
		if (singleInputGridFile == null) {
			singleInputGridFile = new SingleInputGridFile();
			singleInputGridFile.setHistoryKey("PHYML_INPUT_FILE");
			singleInputGridFile.setBorder(new TitledBorder(new LineBorder(new Color(184, 207, 229)), "phyML input file", TitledBorder.LEADING, TitledBorder.TOP, null, new Color(51, 51, 51)));
			if(file!=null){
				try {
					singleInputGridFile.setInputFileUrl(file);
				} catch (RemoteFileSystemException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		return singleInputGridFile;
	}
	private Cpus getCpus() {
		if (cpus == null) {
			cpus = new Cpus();
			cpus.setHistoryKey("PHYML_CPUS");
		}
		return cpus;
	}
	
	public int getCpuValue() {
		return getCpus().getCpus();
	}
	
	public GridFile getPhyMLInputFile() {
		return getSingleInputGridFile().getInputFile();
	}
	
	public String getParameters() {
		return getTextCombo().getText();
	}
	
	public int getWalltimeInSeconds() {
		return getWalltime().getWalltimeInSeconds();
	}
	
	private Walltime getWalltime() {
		if (walltime == null) {
			walltime = new Walltime();
			walltime.setHistoryKey("PHYML_WALLTIME");
		}
		return walltime;
	}

	public void setServiceInterface(ServiceInterface si) {
		// needs to be called so the remote filesystems can be accessed
		getSingleInputGridFile().setServiceInterface(si);
	}

	public void setCommand(String command) {
		this.command=command;
	}

	public void setFile(String file) {
		this.file=file;
	}
}
