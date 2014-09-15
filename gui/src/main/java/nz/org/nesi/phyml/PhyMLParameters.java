package nz.org.nesi.phyml;

import grisu.frontend.view.cli.GrisuCliParameters;

import com.beust.jcommander.Parameter;

public class PhyMLParameters extends GrisuCliParameters {

	@Parameter(names = { "-f", "--input-file" }, description = "the path to the phyml input file")
	private String file;

	@Parameter(names = { "-p", "--phyml-parameters"}, description = "additional phyml parameters")
	private String phyMLParameters;
	
	@Parameter(names = { "-w", "--waltime"}, description = "walltime in seconds or format: 2d10h40m" )
	private String walltime = "60";
	
	@Parameter(names = { "-c", "--cpus"}, description = "how many cpus to use")
	private int cpus = 3;
	
	public String getFile() {
		return file;
	}
	
	public String getPhyMLParameters() {
		return phyMLParameters;
	}
	
	public String getWalltime() {
		return walltime;
	}
	
	public int getCpus() {
		return cpus;
	}
}
