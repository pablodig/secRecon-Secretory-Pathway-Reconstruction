class iCHOSEC_Builder:
    def __init__(self, PSI, protName, sequence, L):
        self.PSI = PSI
        self.protName = protName
        self.sequence = sequence
        self.L = L
        self.rxns = []
        self.rxnNames = []
        self.GPRs = []
        self.connector = None
        self.location = None
        self.clathrin_coeff = None


    def count_AAs(self, sequence):
        """
        Count the number of each amino acid in the given sequence.

        Parameters:
        sequence (str): The amino acid sequence.

        Returns:
        dict: A dictionary with amino acid counts.
        """

        AAs = ['G', 'A', 'V', 'L', 'I', 'M', 'W', 'F', 'P', 'S', 'T', 'C', 'Y', 'N', 'Q', 'E', 'D', 'K', 'R', 'H']
        AAcounts = {aa: sequence.count(aa) for aa in AAs}
        return AAcounts

    def substitute_AAs_count(self, formula, AAcounts):
        """
        Substitute amino acid counts into a reaction formula template.

        Parameters:
        formula (str): The reaction formula with placeholders for amino acids.
        AAcounts (dict): A dictionary with amino acid counts.

        Returns:
        str: The updated reaction formula with amino acid counts substituted.
        """
        template = "? gly_c + ? ala_L_c + ? val_L_c + ? leu_L_c + ? ile_L_c + ? met_L_c + ? trp_L_c + ? phe_L_c + ? pro_L_c + ? ser_L_c + ? thr_L_c + ? cys_L_c + ? tyr_L_c + ? asn_L_c + ? gln_L_c + ? glu_L_c + ? asp_L_c + ? lys_L_c + ? arg_L_c + ? his_L_c"
        new = ''
        aa_list = ['gly_c', 'ala_L_c', 'val_L_c', 'leu_L_c', 'ile_L_c', 'met_L_c', 'trp_L_c', 'phe_L_c', 'pro_L_c',
                   'ser_L_c', 'thr_L_c', 'cys_L_c', 'tyr_L_c', 'asn_L_c', 'gln_L_c', 'glu_L_c', 'asp_L_c', 'lys_L_c',
                   'arg_L_c', 'his_L_c']
        index = 0
        for character in template:
            if character == "?":
                new += str(AAcounts[aa_list[index]])
                index += 1
            else:
                new += character
        newFormula = formula.replace(template, new)
        return newFormula

    def translate_protein(self, entryID, PSIM, PSIM_entries):
        """
        Generate the translation reaction of a protein given its UniProt ID.

        Parameters:
        entryID (str): The UniProt ID of the protein.
        PSIM (list): List of protein synthesis information.
        PSIM_entries (list): List of UniProt IDs corresponding to PSIM.

        Returns:
        str: The translation reaction formula.
        """
        # Obtain protein sequence and length
        PSI_row = PSIM[PSIM_entries.index(str(entryID))] 
        PSI_row = PSI_row.split('\t')
        sequence = PSI_row[11]
        AAcounts = self.count_AAs(sequence)
        N = len(sequence)  # Number to replace atp's, gtp's, ppi's, pi's, h2o's, amp's, and gdp's
        templateFormula = "? h2o_c + ? atp_c + ? gtp_c + ? gly_c + ? ala_L_c + ? val_L_c + ? leu_L_c + ? ile_L_c + ? met_L_c + ? trp_L_c + ? phe_L_c + ? pro_L_c + ? ser_L_c + ? thr_L_c + ? cys_L_c + ? tyr_L_c + ? asn_L_c + ? gln_L_c + ? glu_L_c + ? asp_L_c + ? lys_L_c + ? arg_L_c + ? his_L_c --> ? h_c + ? amp_c + adp_c + ? pi_c + ? gdp_c + ? ppi_c + XXX_c"
        translationFormula = self.substitute_AAs_count(templateFormula, AAcounts)
        translationFormula = translationFormula.replace("? h2o_c", str(2*N-1) + " h2o_c")
        translationFormula = translationFormula.replace("? h_c", str(2*N-1) + " h_c")
        translationFormula = translationFormula.replace("? ppi_c", str(N) + " ppi_c")
        translationFormula = translationFormula.replace("? pi_c", str(2*N-1) + " pi_c")
        translationFormula = translationFormula.replace("? gdp_c", str(2*N-2) + " gdp_c")
        translationFormula = translationFormula.replace("? amp_c", str(N) + " amp_c")
        translationFormula = translationFormula.replace("? gtp_c", str(2*N-2) + " gtp_c")
        translationFormula = translationFormula.replace("? atp_c", str(N+1) + " atp_c")
        translationFormula = translationFormula.replace("XXX_c", str(entryID) + "_c")
        return translationFormula

    def insert_prot_name_in_rxnFormula(self, formula, entryID):
        """
        Replace the XXX template name with the UniProt ID in the reaction formula.

        Parameters:
        formula (str): The reaction formula with XXX as a placeholder.
        entryID (str): The UniProt ID of the protein.

        Returns:
        str: The updated reaction formula with the UniProt ID.
        """
        newFormula = formula.replace("XXX", str(entryID))
        return newFormula

    def insert_prot_name_in_rxnName(self, rxnAbbrev, entryID):
        """
        Replace the XXX template name with the Reaction abbreviation for the UniProt ID.

        Parameters:
        rxnAbbrev (str): The reaction abbreviation with XXX as a placeholder.
        entryID (str): The UniProt ID of the protein.

        Returns:
        str: The updated reaction abbreviation with the UniProt ID.
        """
        newAbbreviation = str(entryID) + "_" + rxnAbbrev
        return newAbbreviation

    def addPathway(self, pathwayName, listOfRxns, listOfRxnsNames):
        """
        Add reactions and reaction names for a given pathway.

        Parameters:
        pathwayName (str): The name of the pathway.
        listOfRxns (list): The list of reactions to add to.
        listOfRxnsNames (list): The list of reaction names to add to.

        Returns:
        tuple: Updated lists of reactions and reaction names.
        """
        newList = listOfRxns
        newList2 = listOfRxnsNames
        for i in range(len(rxnPathway)):
            if rxnPathway[i] == pathwayName:
                newList.append(rxnFormula[i])
                newList2.append(rxnAbbreviation[i])
        return newList, newList2

    def build_reactions(self):
        PSI = self.PSI
        protName = self.protName
        sequence = self.sequence
        L = self.L
        rxns = self.rxns
        rxnNames = self.rxnNames
        GPRs = self.GPRs

        if PSI[7] == '[x]' or PSI[7] == '[l]' or PSI[7] == '[d]':
            rxns.append(self.connector + ' --> XXX-preClathrin[g]')
            rxnNames.append('Start_Clathrin_vesicle')

            clath_rxns, clath_names = self.add_pathway('Clathrin vesicles', [], [])
            for i in range(len(clath_rxns)):
                clath_rxns[i] = clath_rxns[i].replace("!", str(self.clathrin_coeff))
                rxns.append(clath_rxns[i])
                rxnNames.append(clath_names[i])

            self.connector = 'XXX_mature[cv]'
            self.location = PSI[7]
            rxns.append(self.connector + ' --> XXX_mature' + self.location)
            rxnNames.append('Final_location_' + self.location)
            rxns.append('XXX_mature' + self.location + ' --> ')
            rxnNames.append(protName + '_Final_demand')

        else:
            self.location = PSI[7]
            rxns.append(self.connector + ' --> XXX-preSV[g]')
            rxnNames.append('Start_Secretion')

            sv_rxns, sv_names = self.add_pathway('SV', [], [])
            for i in range(len(sv_rxns)):
                sv_rxns[i] = sv_rxns[i].replace("!", str(self.clathrin_coeff))
                rxns.append(sv_rxns[i])
                rxnNames.append(sv_names[i])

            if self.location == '':
                self.location = '[e]'
            rxns.append('XXX_mature[sv]' + ' --> XXX_mature' + self.location)
            rxnNames.append('Final_location_' + self.location)
            rxns.append('XXX_mature' + self.location + ' --> ')
            rxnNames.append('Final_demand')

        if 'SP_degradation' in rxnNames:
            SPaas = self.count_AAs(sequence[0:22])
            rxns[rxnNames.index('SP_degradation')] = self.substitute_AAs_count(rxns[rxnNames.index('SP_degradation')], SPaas)
        if 'Ubiquitination_degradation' in rxnNames:
            new_aas = self.count_AAs(sequence[22:])
            rxns[rxnNames.index('Ubiquitination_degradation')] = rxns[rxnNames.index('Ubiquitination_degradation')].replace("!", "?")
            rxns[rxnNames.index('Ubiquitination_degradation')] = self.substitute_AAs_count(rxns[rxnNames.index('Ubiquitination_degradation')], new_aas)
            rxns[rxnNames.index('Ubiquitination_degradation')] = rxns[rxnNames.index('Ubiquitination_degradation')].replace("?", str(L-22))

        GPRs = self.getGPRsFromRxnNames(rxnNames)

        for i in range(len(rxns)):
            rxns[i] = self.insert_prot_name_in_rxnFormula(rxns[i], PSI[0])
        for i in range(len(rxnNames)):
            rxnNames[i] = self.insert_prot_name_in_rxnName(rxnNames[i], PSI[0])

        rxns, rxnNames, GPRs = self.addCanonicalRxns(rxns, rxnNames, GPRs)
        return rxns, rxnNames, GPRs
