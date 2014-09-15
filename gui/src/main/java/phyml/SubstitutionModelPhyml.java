package phyml;

import javax.swing.*;
import java.awt.*;

/**
 * JPanel implementing all components necessary for specifying the
 * "Substitution Model" parameters.
 *
 * @author Christoph Knapp
 */

public class SubstitutionModelPhyml extends JPanel {
	/**
	 * default id
	 */
	private static final long serialVersionUID = 1L;
	static SubModel sM;
	private OptimiseFreqAndRatio oFAR;
	private ProportionAndSubrate pAS;
	private NumSubCatGammaAverage nSCGA;

	/**
	 * JPanel implementing all components necessary for specifying the
	 * "Substitution Model" parameters.
	 *
	 * @param molecularType
	 *            : whether DNA or AA is "Data Type"
	 */
	public SubstitutionModelPhyml(String molecularType) {
		CustomGridLayout layout = new CustomGridLayout();
		setLayout(layout);
		layout.setDimensions(1, 0.1);
		add(new Separator("Substitution Model", Color.lightGray, null,false));
		layout.setDimensions(1, 0.4);
		sM = new SubModel(molecularType);
		add(sM);
		layout.setDimensions(1, 0.08);
		oFAR = new OptimiseFreqAndRatio("DNA");
		add(oFAR);
        sM.setOptimiseFreqAndRatioPanel(oFAR);
		layout.setDimensions(1, 0.16);
		pAS = new ProportionAndSubrate();
		add(pAS);
        sM.setProportionAndSubratePanel(pAS);
		layout.setDimensions(1, 0.24);
		nSCGA = new NumSubCatGammaAverage();
		add(nSCGA);
        sM.setNumSubCatGammaAveragePanel(nSCGA);
	}

	/**
	 * Retrieves the type of Substitution model from a SubModel object.
	 *
	 * @return -m (or --model) model name<br>
	 *         model name : substitution model name.<br>
	 *         - Nucleotide-based models: HKY85 (default) | JC69 | K80 | F81 |
	 *         F84 | TN93 | GTR | custom. The custom option can be used to define
	 *         a new substitution model. A string of six digits identifies the
	 *         model. For instance, 000000 corresponds to F81 (or JC69 provided
	 *         the distribution of nucleotide frequencies is uni-form). 012345
	 *         corresponds to GTR. This option can be used for encoding any
	 *         model that is a nested within GTR. See Section 6.1.2. NOTE: the
	 *         substitution parameters of the custom model will be optimised so
	 *         as to maximise the likelihood. It is possible to specify and fix
	 *         (i.e., avoid optimisation) the values of the substitution rates
	 *         only through the PHYLIP-like interface.<br>
	 *         - Amino-acid based models: LG (default) WAG | JTT | MtREV |
	 *         Dayhoff | DCMut | RtREV | CpREV | VT | Blosum62 | MtMam | MtArt |
	 *         HIVw | HIVb | custom The custom option is useful when one wants
	 *         to use an amino-acid substitution model that is not available by
	 *         default in PhyML. The symmetric part of the rate matrix, as well
	 *         as the equilibrium amino-acid frequencies, are given in a file
	 *         which name is asked for by the program. The format of this file
	 *         is described in section 7.4 in the phyml manual.
	 */
	public String getSubModel() {
		return sM.getSubModel();
	}

	/**
	 * Retrieves the path to the Amino Acid rate file specified by the user.
	 *
	 * @return --aa rate file file name This option is compulsory when analysing
	 *         amino-acid sequences under a ‘custom’ model. file name should
	 *         provide a rate matrix and equilibrium amino acid in PAML format
	 *         (see Section ).
	 */
	public String getAARateFileName() {
		return sM.getAARateFileName();
	}

	/**
	 * Retrieve the frequencies if custom model was selected and nucleotide
	 * frequencies were specified.
	 *
	 * @return -f e, m, or “fA,fC,fG,fT”<br>
	 *         Nucleotide or amino-acid frequencies.<br>
	 *         - e : the character frequencies are determined as follows : -
	 *         Nucleotide sequences: (Empirical) the equilibrium base
	 *         frequencies are estimated by counting the occurence of the
	 *         different bases in the alignment. - Amino-acid sequences:
	 *         (Empirical) the equilibrium amino-acid frequencies are estimated
	 *         by counting the occurence of the different amino-acids in the
	 *         alignment.<br>
	 *         - m : the character frequencies are determined as follows : -
	 *         Nucleotide sequences: (ML) the equilibrium base frequencies are
	 *         estimated using maximum likelihood. - Amino-acid sequences:
	 *         (Model) the equilibrium amino-acid frequencies are estimated
	 *         using the frequencies defined by the substitution model.<br>
	 *         - “fA,fC,fG,fT” : only valid for nucleotide-based models. fA, fC,
	 *         fG and fT are floating numbers that correspond to the frequencies
	 *         of A, C, G and T respectively.<br>
	 */
	public String getFrequencies() {
		return sM.getFrequencies();
	}

	/**
	 * Retrieves the TS/TV ratio for the given Substitution Model.
	 *
	 * @return -t (or --ts/tv) ts/tv ratio ts/tv ratio: transition/transversion
	 *         ratio. DNA sequences only. Can be a fixed positive value (e.g.,
	 *         4.0) or type e to get the maximum likelihood estimate.
	 */
	public String getTsTvRatio() {
		return oFAR.getTsTvRatio();
	}

	/**
	 * Retrieves the proportion of invariable sites parameter from the
	 * ProportionAndSubrate panel.
	 *
	 * @return -v (or --pinv) prop invar prop invar: proportion of invariable
	 *         sites. Can be a fixed value in the [0,1] range or type e to get
	 *         the maximum likelihood estimate.
	 */
	public String getPropInvarSites() {
		return pAS.getPropInvarSites();
	}

	/**
	 * Retrieves the number of substitution categories from
	 * NumSubCatGammaAverage object.
	 *
	 * @return -c (or --nclasses) nb subst cat nb subst cat: number of relative
	 *         substitution rate categories. Default: nb subst cat=4. Must be a
	 *         positive integer.
	 */
	public String getNumSubCat() {
		return nSCGA.getNumSubCat();
	}

	/**
	 * Retrieves the Alpha value as a String from NumSubCatGammaAverage object.
	 *
	 * @return -a (or --alpha) gamma gamma: value of the gamma shape parameter.
	 *         Can be a fixed positive value or e to get the maximum likelihood
	 *         estimate. The value of this parameter is estimated in the maximum
	 *         likelihood framework by default.
	 */
	public String getAlpha() {
		return nSCGA.getAlpha();
	}

	/**
	 * Retrieves whether the median or the mean is used as average from a
	 * NumSubCatGammaAverage object.
	 *
	 * @return --use median The middle of each substitution rate class in the
	 *         discrete gamma distribution is taken as the median. The mean is
	 *         used by default.
	 */
	public String getUseMedian() {
		return nSCGA.getUseMedian();
	}

	/**
	 * Retrieves the optimised rate parameter from a SubModel object.
	 *
	 * @return --free rates As an alternative to the discrete gamma model, it is
	 *         possible to estimate the (relative) rate in each class of the
	 *         (mixture) model and the corresponding frequencies. This model has
	 *         more parameters than the discrete gamma one but usually provides
	 *         a significantly better fit to the data.
	 */
	public boolean isOptimisedRateParameter() {
		return sM.isOptimisedRateParameter();
	}

	/**
	 * Sets all components in the NumSubCatGammaAverage object to visible or not
	 * visible.
	 *
	 * @param b
	 *            boolean : true if visible, false otherwise.
	 */
	public void setNumSubCatGammaAverageCompVisible(boolean b) {
		nSCGA.setCompVisible(b);
	}

	/**
	 * Forwards the change of the molecule type to an OptimiseFreqAndRatio
	 * object.
	 *
	 * @param string
	 *            String : either "DNA" or "AA"
	 */
	public void setMoleculeType(String string) {
		oFAR.setMoleculeType(string);
		sM.setMoleculeType(string);
	}
	/**
	 * Enables or disables the dropdownmenu for changing the Transition/Transversion ratio.
	 *
	 * @param off
	 * boolean : if true disabled otherwise enabled
	 */
	public void setRatioBoxOff(boolean off) {
		oFAR.setRatioBoxOff(off);
	}

	public void setRatioBoxDefault() {
		oFAR.setRatioBoxDefault();
	}
}
