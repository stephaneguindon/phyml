package nz.org.nesi.phyml.swing;

import com.jgoodies.forms.layout.ColumnSpec;
import com.jgoodies.forms.layout.FormLayout;
import com.jgoodies.forms.layout.FormSpecs;
import com.jgoodies.forms.layout.RowSpec;
import grisu.control.ServiceInterface;
import grisu.frontend.control.jobMonitoring.RunningJobManagerManager;
import grisu.frontend.model.job.GrisuJob;
import grisu.frontend.view.swing.jobcreation.JobCreationPanel;
import grisu.frontend.view.swing.jobcreation.widgets.SubmissionLogPanel;
import grisu.model.dto.GridFile;
import nz.org.nesi.phyml.Client;
import org.apache.commons.lang.StringUtils;
import org.jdesktop.swingx.JXErrorPane;
import org.jdesktop.swingx.error.ErrorInfo;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.logging.Level;

public class PhyMLJobCreationPanel extends JPanel implements JobCreationPanel {
	// the button to submit a job
	private JButton btnSubmit;
	// a ready-made widget that, once connected to a GrisuJob, tracks the
	// progress of a job submission...
	private SubmissionLogPanel submissionLogPanel;

	// the serviceinterface. gets set when the setServiceInterface method is
	// called. we need it in order to submit the job and do everything else...
	private ServiceInterface si;
	private PhyMLParameterInputMask phyMLParameterInputMask;
	private String command;
	private String file;
	public PhyMLJobCreationPanel() {
		setLayout(new FormLayout(new ColumnSpec[] {
				FormSpecs.RELATED_GAP_COLSPEC,
				ColumnSpec.decode("max(112dlu;default):grow"),
				FormSpecs.RELATED_GAP_COLSPEC,},
			new RowSpec[] {
				FormSpecs.DEFAULT_ROWSPEC,
				FormSpecs.RELATED_GAP_ROWSPEC,
				FormSpecs.DEFAULT_ROWSPEC,
				FormSpecs.RELATED_GAP_ROWSPEC,
				RowSpec.decode("default:grow"),
				FormSpecs.RELATED_GAP_ROWSPEC,}));
		add(getPhyMLParameterInputMask(), "2, 1");
		add(getBtnSubmit(), "2, 3, right, center");
		add(getSubmissionLogPanel(), "2, 5, fill, fill");
	}

	@Override
	public boolean createsBatchJob() {
		// we are only creating a single job with this panel
		return false;
	}

	@Override
	public boolean createsSingleJob() {
		// yep
		return true;
	}

	// creating the submit button and connecting it with the submit action
	private JButton getBtnSubmit() {
		if (btnSubmit == null) {
			btnSubmit = new JButton("Submit");
			btnSubmit.setAlignmentX(Component.RIGHT_ALIGNMENT);
			btnSubmit.addActionListener(new ActionListener() {
				@Override
				public void actionPerformed(ActionEvent arg0) {
					submitJob();
				}
			});
		}
		return btnSubmit;
	}

	// usually you just return the object itself
	@Override
	public JPanel getPanel() {
		return this;
	}

	// this is what is going to be displayed in the client on the upper left
	// side
	@Override
	public String getPanelName() {
		return "phyML";
	}

	// creating a submission log panel
	private SubmissionLogPanel getSubmissionLogPanel() {
		if (submissionLogPanel == null) {
			submissionLogPanel = new SubmissionLogPanel();
			submissionLogPanel.setAlignmentX(Component.RIGHT_ALIGNMENT);
		}
		return submissionLogPanel;
	}

	// here you return the name of the application package that your client
	// submits jobs for.
	// e.g. Java, Python, MrBayes, ....
	@Override
	public String getSupportedApplication() {
		return "phyml";
	}

	// convenience method to disable and re-enable the button while the
	// submission is ongoing
	private void lockUI(final boolean lock) {

		SwingUtilities.invokeLater(new Thread() {

			@Override
			public void run() {
				getBtnSubmit().setEnabled(!lock);
			}

		});

	}

	// here we collect the reference to the serviceinterface. Also, usually
	// here's where you do some
	// initialization of your panel (if you need to). E.g. populating a combobox
	// with a list of submission locations or versions of an application
	// package...
	@Override
	public void setServiceInterface(ServiceInterface si) {
		this.si = si;
		getPhyMLParameterInputMask().setServiceInterface(si);
	}

	// the method to actually submit the job
	private void submitJob() {

		// we submit the job in it's own thread in order to not block the swing
		// ui
		new Thread() {
			@Override
			public void run() {
				try {
					// disable the submit button so the user can't inadvertently
					// submit 2 jobs in a row
					lockUI(true);

					// get the values the user selected
					String parameters = getPhyMLParameterInputMask().getParameters();
					int cpus = getPhyMLParameterInputMask().getCpuValue();
					int walltime = getPhyMLParameterInputMask().getWalltimeInSeconds();

					GridFile inputFile = null;
					try {
						inputFile = getPhyMLParameterInputMask().getPhyMLInputFile();
					} catch (Exception e) {
						throw new Exception("No valid input file specified.");
					}

					// do some validation
					if ( inputFile == null ) {
						throw new Exception("No input file specified.");
					}
					if ( StringUtils.isBlank(parameters) ) {
						parameters = Client.DEFAULT_PHYML_PARAMETERS;
					}


					// now, let's create the job
					GrisuJob job = new GrisuJob(si);
					// ... and connect it to the submission log panel so the
					// user can see that there's something going on...
					getSubmissionLogPanel().setJobObject(job);
					// .. every job needs its own unique jobname. In this case we'll use the name of the input file and a timestamp

					job.setApplication("phyml");
					job.setApplicationVersion("20120412");
					job.setSubmissionLocation("pan:pan.nesi.org.nz");


					job.setTimestampJobname(inputFile.getName());
					// this is required and the most important part of the job
					// creation process. Grisu needs to know what job you want
					// to submit, after all...
					job.setCommandline("/share/apps/phyml-20120412/bin/phyml-mpi -i "+inputFile.getName()+" "+parameters);
					// setting the walltime of a job is also important
					job.setWalltimeInSeconds(walltime);
					// setting the cpus
					job.setCpus(cpus);

					// finally, we also have to attach the input file so it gets uploaded / copied into the job directory
					job.addInputFileUrl(inputFile.getUrl());

					// using the '/nz/nesi' group hard-coded here, we might need to change that later
					// creating the job using the RunningJobManager notifies the job monitoring panel that there is a new job
					RunningJobManagerManager.getDefault(si).createJob(job, "/nz/nesi");

					// last not least, we stage in files and submit the job
					job.submitJob();

				} catch (Exception e) {
					// if something goes wrong, we want to show the user
					ErrorInfo info = new ErrorInfo("Job submission error",
							"Can't submit job:\n\n" + e.getLocalizedMessage(),
							null, "Error", e, Level.SEVERE, null);

					JXErrorPane pane = new JXErrorPane();
					pane.setErrorInfo(info);
					// the following line could be used to show a button to
					// submit
					// a bug/ticket to a bug tracking system
					// pane.setErrorReporter(new GrisuErrorReporter());

					JXErrorPane.showDialog(
							PhyMLJobCreationPanel.this.getRootPane(), pane);
				} finally {
					// enable the submit button again
					lockUI(false);
				}

			}
		}.start();
	}
	private PhyMLParameterInputMask getPhyMLParameterInputMask() {
		if (phyMLParameterInputMask == null) {
			phyMLParameterInputMask = new PhyMLParameterInputMask();
			phyMLParameterInputMask.setCommand(command);
			phyMLParameterInputMask.setFile(file);
		}
		return phyMLParameterInputMask;
	}

	public void setCommand(String command) {
		this.command=command;
	}

	public void setFile(String file) {
		this.file=file;
	}
}
