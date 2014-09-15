package nz.org.nesi.phyml;

import grisu.control.ServiceInterface;
import grisu.control.exceptions.JobPropertiesException;
import grisu.control.exceptions.JobSubmissionException;
import grisu.frontend.control.jobMonitoring.RunningJobManagerManager;
import grisu.frontend.model.job.GrisuJob;
import grisu.frontend.view.cli.GrisuCliClient;
import grisu.jcommons.constants.Constants;
import grisu.jcommons.utils.WalltimeUtils;
import grisu.model.FileManager;
import org.apache.commons.lang.StringUtils;
import phyml.StandardOutPanel;

public class Client extends GrisuCliClient<PhyMLParameters> {

	public static final String DEFAULT_PHYML_PARAMETERS = "-d nt -q -n 1 -m HKY85 -f e -t e -v 0.0 -c 4 -a e -o tl --print_site_lnl --print_trace --no_memory_check";

	public Client(PhyMLParameters params, String[] args) throws Exception {
		super(params, args);
	}

	@Override
	public void run() {

		String phyml_input_file_path = getCliParameters().getFile();
		if ( StringUtils.isBlank(phyml_input_file_path) ) {
			StandardOutPanel.setInput("No input file specified. Exiting...");
			//System.exit(1);
		}
		String phyml_input_filename = FileManager.getFilename(phyml_input_file_path);
		StandardOutPanel.setInput("File to use for the job: " + phyml_input_file_path);

		String parameters = getCliParameters().getPhyMLParameters();
		StandardOutPanel.setInput("Using phymlparameters: "+parameters);
		if ( StringUtils.isBlank(parameters) ) {
			StandardOutPanel.setInput("Using default phymlparameters: "+DEFAULT_PHYML_PARAMETERS);
		}

		String walltimeTemp = getCliParameters().getWalltime();
		int walltime = -1;
		try {
			walltime = WalltimeUtils.fromShortStringToSeconds(walltimeTemp);
			StandardOutPanel.setInput("Walltime in seconds: "+walltime);
		} catch (Exception e1) {
			StandardOutPanel.setInput("Invalid walltime specified: "+walltimeTemp);
			System.exit(1);
		}

		int cpus = getCliParameters().getCpus();
		StandardOutPanel.setInput("CPUs to use: : "+cpus);

		// ============ Login =========================================================
		// all login stuff is implemented in the parent class
		StandardOutPanel.setInput("Getting serviceinterface...");
		ServiceInterface si = null;
		try {
			si = getServiceInterface();
		} catch (Exception e) {
			StandardOutPanel.setInput("Could not login: " + e.getLocalizedMessage());
			System.exit(1);
		}


		// ============ Job creation ======================================================
		StandardOutPanel.setInput("Creating job...");
		GrisuJob job = new GrisuJob(si);

		// we can savely hard-code that. might need a config option for version though
		job.setApplication("phyml");
		job.setApplicationVersion("20120412");
		job.setSubmissionLocation("pan:pan.nesi.org.nz");

		job.setCommandline("/share/apps/phyml-20120412/bin/phyml-mpi -i "+phyml_input_filename+" "+parameters);
		job.addInputFileUrl(phyml_input_file_path);

		job.setCpus(cpus);
		job.setWalltimeInSeconds(walltime);

		job.setTimestampJobname(phyml_input_filename);

		StandardOutPanel.setInput("Set jobname to be: " + job.getJobname());

		try {
			StandardOutPanel.setInput("Creating job on backend...");
			// group might need to change
			RunningJobManagerManager.getDefault(si).createJob(job, "/nz/nesi");
		} catch (JobPropertiesException e) {
			StandardOutPanel.setInput("Could not create job: "
					+ e.getLocalizedMessage());
			System.exit(1);
		}

		// ============ Job submission ======================================================
		try {
			StandardOutPanel.setInput("Submitting job to the grid...");
			job.submitJob();
		} catch (JobSubmissionException e) {
			StandardOutPanel.setInput("Could not submit job: "
					+ e.getLocalizedMessage());
			System.exit(1);
		} catch (InterruptedException e) {
			StandardOutPanel.setInput("Jobsubmission interrupted: "
					+ e.getLocalizedMessage());
			System.exit(1);
		}

		StandardOutPanel.setInput("Job submission finished.");
		StandardOutPanel.setInput("Job submitted to: "
				+ job.getJobProperty(Constants.SUBMISSION_SITE_KEY));


		// ============ Waiting for job to finish ================================================
		StandardOutPanel.setInput("Waiting for job to finish...");

		// for a realy workflow, don't check every 5 seconds since that would
		// put too much load on the backend/gateways
		job.waitForJobToFinish(5);

		// ============ Job results ==============================================================
		StandardOutPanel.setInput("Job finished with status: "
				+ job.getStatusString(false));

		StandardOutPanel.setInput("Stdout: " + job.getStdOutContent());
		try {
			StandardOutPanel.setInput("Stderr: " + job.getStdErrContent());
		} catch (Exception e){
			// no stderr file, no worries...
		}

		//System.exit(0);
	}

}
