package nz.org.nesi.phyml.swing;

import grisu.control.ServiceInterface;
import grisu.frontend.control.jobMonitoring.RunningJobManagerManager;
import grisu.frontend.model.job.GrisuJob;
import grisu.model.FileManager;
import phyml.PhymlPanel;
import phyml.StandardOutPanel;

import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.List;

public class PhyMLJobSubmit {

	private final ServiceInterface si;
	private final String command;
	private final String file;

	public PhyMLJobSubmit(ServiceInterface si, String command, String file) {
		this.si = si;
		this.command = command;
		this.file = file;
	}

	public void submit() throws Exception {

		GrisuJob job = new GrisuJob(si);


		job.setApplication("PhyML");
		job.setApplicationVersion("20120412-sandybridge");
		job.setSubmissionLocation("pan:gram.uoa.nesi.org.nz");

		//job.setCommandline("phyml-mpi "+command);
        job.setForce_mpi(true);

        job.setHostCount(0);
        job.setCpus(10);

        job.setCommandline("phyml-mpi "+command);
		job.addInputFileUrl(file);

        job.addEnvironmentVariable("LL_VAR", "requirements=(Feature==\"sandybridge\")");
        job.addEnvironmentVariable("MPI_PREFIX", "mpiexec.hydra --machinefile \\${LOADL_HOSTFILE} -genv I_MPI_FABRICS dapl -genv I_MPI_PIN_PROCESSOR_LIST=\"grain=cache2,shift=sock\" -envall");

		try {
			job.setWalltime("10m");
		} catch (Exception e) {
			throw new RuntimeException(e);
		}

		job.setTimestampJobname(FileManager.getFilename(file));

		PropertyChangeListener l = new PropertyChangeListener() {

			@Override
			public void propertyChange(PropertyChangeEvent evt) {

				String propName = evt.getPropertyName();
				if ("submissionLog".equals(propName)) {
					final List<String> log = (List<String>) evt.getNewValue();
					String text = log.get(log.size() - 1);
					StandardOutPanel.setInput(text);
				}

			}
		};

		job.addPropertyChangeListener(l);

		try {
		RunningJobManagerManager.getDefault(si).createJob(job, "/nz/nesi");

		job.submitJob();

		} finally {
			job.removePropertyChangeListener(l);
		}


        PhymlPanel.SetSubmit(true);
        PhymlPanel.loadTrees();

	}

}
