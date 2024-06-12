class iCHOSEC_Builder:
    def __init__(self, entryID, PSIM, rxnPathway, rxnFormula, rxnAbbreviation, rxnConditions, rxnGPR, rxnComponents):
        """
        Initialize the iCHOSEC_Builder class with the provided parameters.

        Parameters:
        entryID (int or str): The ID of the protein entry.
        PSIM (DataFrame): Protein Synthesis Information matrix.
        rxnPathway (str): The reaction pathway.
        rxnFormula (str): The reaction formula.
        rxnAbbreviation (str): The reaction abbreviation.
        rxnConditions (str): The reaction conditions.
        rxnGPR (str): The gene-protein-reaction association.
        """
        self.entryID = entryID
        self._validate_input(PSIM, entryID)

        self.PSI = PSIM.loc[entryID]
        self.protName = self.PSI.iloc[0]
        self.sequence = self.PSI.iloc[10]
        self.L = float(self.PSI.iloc[1]) # Protein Length
        self.MW = float(self.PSI.iloc[2])  # Molecular weight
        self.SP = str(self.PSI.iloc[3])  # Signal Peptide
        self.DSB = str(self.PSI.iloc[4])  # Disulphide Bond
        self.GPI = str(self.PSI.iloc[5])  # GPI
        self.NG = str(self.PSI.iloc[6])  # N-glycosylation
        self.OG = str(self.PSI.iloc[7])  # O-glycosylation
        self.TMD = str(self.PSI.iloc[8])  # Transmembrane Domains
        self.subloc = str(self.PSI.iloc[9])  # Subcellular Localization
        self.rxns = []
        self.rxnNames = []
        self.GPRs = []
        self.connector = None
        self.location = None
        self.clathrin_coeff = None
        self.rxnPathway = rxnPathway
        self.rxnFormula = rxnFormula
        self.rxnAbbreviation = rxnAbbreviation
        self.rxnConditions = rxnConditions
        self.rxnGPR = rxnGPR
        self.rxnComponents = rxnComponents

        print(f'''
            Protein: {self.protName}
            Length: {self.L}
            MW: {self.MW}
            Signal Peptide: {self.SP}
            Disulfide Bonds: {self.DSB}
            ''')

    def _validate_input(self, PSIM, entryID):
        """
        Validate the input parameters.

        Parameters:
        PSIM (DataFrame): Protein Synthesis Information matrix.
        entryID (int or str): The ID of the protein entry.

        Raises:
        ValueError: If PSIM does not contain the entryID or required columns.
        """
        if entryID not in PSIM.index:
            raise ValueError(f"entryID {entryID} not found in PSIM")

        required_columns = [
            'Protein names', 'Length', 'Mass', 'SP', 'DSB', 'GPI', 
            'NG', 'OG', 'TMD', 'Location', 'Sequence'
        ]
        for column in required_columns:
            if column not in PSIM.columns:
                raise ValueError(f"PSIM is missing required column: {column}")

    def count_AAs(self, sequence):
        """
        Count the number of each amino acid in the given sequence.

        Parameters:
        sequence (str): The amino acid sequence.

        Returns:
        list: A list with amino acid counts.
        """
        AAcounts = []
        AAs = ['G', 'A', 'V', 'L', 'I', 'M', 'W', 'F', 'P', 'S', 'T', 'C', 'Y', 'N', 'Q', 'E', 'D', 'K', 'R', 'H']
        for aa in AAs:
            AAcounts.append(sequence.count(aa))
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
        
        index = 0
        for character in template:
            if character == "?":
                new += str(AAcounts[index])
                index += 1
            else:
                new += character
        newFormula = formula.replace(template, new)
        return newFormula

    def translate_protein(self):
        """
        Generate the translation reaction of a protein given its UniProt ID.

        Returns:
        str: The translation reaction formula.
        """
        
        AAcounts = self.count_AAs(self.sequence)
        N = len(self.sequence)  # Number to replace atp's, gtp's, ppi's, pi's, h2o's, amp's, and gdp's
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
        translationFormula = translationFormula.replace("XXX_c", str(self.entryID) + "_c")
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
        if rxnAbbrev.startswith("SK_"): # exclude Sink reactions for Sec Pathway components
            processed_reaction = rxnAbbrev
        else:
            processed_reaction = str(entryID) + "_" + rxnAbbrev
        return processed_reaction

    def addPathway(self, pathwayName, listOfRxns, listOfRxnsNames, listOfComponents):
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
        for i in range(len(self.rxnPathway)):
            if self.rxnPathway[i] == pathwayName:
                newList.append(self.rxnFormula[i])
                newList2.append(self.rxnAbbreviation[i])
                # Add Sink reactions for the components of each reaction added
                if listOfComponents[i].split(",") != ['']:
                    for component in listOfComponents[i].split(","):
                        newList.append(component + ' <=>')
                        newList2.append('SK_' + component)

        return newList, newList2

    def getGPRsFromRxnNames(self, listOfRxnsNames):
        """
        Create a list of GPRs given a list of reaction names.

        Parameters:
        listOfRxnsNames (list): List of reaction names.

        Returns:
        list: List of GPR associations.
        """
        GPR_list = []
        for reaction in listOfRxnsNames:
            if reaction in self.rxnAbbreviation:
                GPR_list.append(self.rxnGPR[self.rxnAbbreviation.index(reaction)])
                if self.rxnGPR[self.rxnAbbreviation.index(reaction)] == None:
                    print(reaction, self.rxnAbbreviation.index(reaction), self.rxnGPR[self.rxnAbbreviation.index(reaction)])
            else:
                GPR_list.append('')
        return GPR_list

    def addCanonicalRxns(self, listOfRxns, listOfRxnNames, listOfGPRs):
        """
        Add canonical reactions to the overall list of reactions, reaction names, and GPRs.

        Parameters:
        listOfRxns (list): The list of reactions.
        listOfRxnNames (list): The list of reaction names.
        listOfGPRs (list): The list of GPR associations.

        Returns:
        tuple: Updated lists of reactions, reaction names, and GPRs.
        """
        newListRxns = listOfRxns
        newListNames = listOfRxnNames
        newListGPRs = listOfGPRs

        # Add canonical reactions
        newListRxns, newListNames = self.addPathway("Canonical", listOfRxns, listOfRxnNames, self.rxnComponents)

        # Add canonical GPRs
        r = []
        n = []
        r, n = self.addPathway("Canonical", r, n, self.rxnComponents)
        newGPRs = self.getGPRsFromRxnNames(n)
        for gpr in newGPRs:
            newListGPRs.append(gpr)

        return newListRxns, newListNames, newListGPRs

    def addPathwayFromCondition(self, conditionName, listOfRxns, listOfRxnsNames):
        """
        Add reactions and reaction names for a given condition.
    
        Parameters:
        conditionName (str): The condition name.
        listOfRxns (list): The list of reactions to add to.
        listOfRxnsNames (list): The list of reaction names to add to.
    
        Returns:
        tuple: Updated lists of reactions and reaction names.
        """
        newList = listOfRxns
        newList2 = listOfRxnsNames
        for i in range(len(self.rxnConditions)):
            if self.rxnConditions[i] == conditionName:
                newList.append(self.rxnFormula[i])
                newList2.append(self.rxnAbbreviation[i])
        return newList, newList2


    def generateProteinSpecificRxns_A(self):
        """
        Generate protein-specific reactions list based on the entry ID already defined in the class.

        Returns:
        tuple: Lists of reactions, reaction names, and GPRs.
        """
        
        Kv = 0.7
        V = self.MW * 1.21 / 1000.0  # Protein Volume in nm^3
        self.clathrin_coeff = int(round(29880.01 * Kv / V))  # Number of proteins per clathrin vesicle
        copi_coeff = int(round(143793.19 * Kv / V))
        copii_coeff = int(round(268082.35 * Kv / V))
        self.connector = ''

        # Prepare vectors that will store reactions and components
        self.rxns = []
        self.rxnNames = []

        # Add translation reaction
        translation_reaction = self.translate_protein()
        self.rxns.append(translation_reaction)
        self.rxnNames.append("TRANSLATION_protein")

        # Add translocation reactions SP = "Signal Peptide"
        if self.SP == '0':  # If it doesn't have signal peptide then ignore
            self.GPRs = self.getGPRsFromRxnNames(self.rxnNames)
            # Change the reaction names and their formulas to include the UniProt Identifier
            for i in range(len(self.rxns)):
                self.rxns[i] = self.insert_prot_name_in_rxnFormula(self.rxns[i], self.entryID)
            for i in range(len(self.rxnNames)):
                self.rxnNames[i] = self.insert_prot_name_in_rxnName(self.rxnNames[i], self.entryID)
            self.rxns.append(self.entryID + '_c --> ')
            self.rxnNames.append(self.entryID + '_Final_demand')
            self.GPRs.append('')
            return self.rxns, self.rxnNames, self.GPRs

        elif self.SP == '1':  # Translocate protein if it has signal peptide
            if self.L <= 160:
                self.rxns, self.rxnNames = self.addPathway("Post-translational Translocation", self.rxns, self.rxnNames)
                if self.subloc == '[e]' or self.subloc == '':
                    self.rxns, self.rxnNames = self.addPathway("Post-translational Translocation (Secretory protein)", self.rxns, self.rxnNames)
                if self.subloc == '[pm]' or self.TMD != '0':
                    self.rxns, self.rxnNames = self.addPathway("Post-translational Translocation (Tail anchored membrane protein)", self.rxns, self.rxnNames)
            else:
                self.rxns, self.rxnNames = self.addPathway("Translocation", self.rxns, self.rxnNames, self.rxnComponents)

            number_BiP = self.L / 40  # Number of BiPs depends on protein length
            for i in range(len(self.rxns)):
                self.rxns[i] = self.rxns[i].replace("!", str(number_BiP))
            self.connector = 'XXX_r'

        # Add Disulphide Bond Reactions
        if self.DSB != '0':
            number_DSB = str(self.DSB)
            DSBrxns = []
            DSBrxnNames = []
            DSBrxns.append(self.connector + ' --> XXX_preDSB_r')
            DSBrxnNames.append('Start_DSB')
            DSBrxns, DSBrxnNames = self.addPathwayFromCondition('DSB>0', DSBrxns, DSBrxnNames)
            for i in range(len(DSBrxns)):
                if '?' in DSBrxns[i]:
                    if number_DSB == '1':
                        DSBrxns[i] = DSBrxns[i].replace('? ', '')
                    else:
                        DSBrxns[i] = DSBrxns[i].replace('?', number_DSB)

            self.connector = 'XXX_DSB_r'

            for reaction in DSBrxns:
                self.rxns.append(reaction)
            for reactionName in DSBrxnNames:
                self.rxnNames.append(reactionName)

        # Add GPI reactions
        if self.GPI == '1':
            self.rxns.append(self.connector + ' --> XXX_preGPI_r')
            self.rxnNames.append('Start_GPI')
            self.rxns, self.rxnNames = self.addPathwayFromCondition('GPI=1', self.rxns, self.rxnNames)
            self.connector = 'XXX-dgpi_cho_r'

        # Add N-glycosylation reactions
        if self.NG != '0':  # Has N-glycans?
            self.rxns.append(self.connector + ' --> XXX_preNG_r')
            self.rxnNames.append('Start_NG')
            number_Nglycans = str(self.NG)  # Get number of N-Glycans
            NGlyrxns = []
            NGlyrxnNames = []

            NGlyrxns, NGlyrxnNames = self.addPathwayFromCondition('NG>0', NGlyrxns, NGlyrxnNames)
            NGlyrxns, NGlyrxnNames = self.addPathway('Golgi processing N', NGlyrxns, NGlyrxnNames)
            for i in range(len(NGlyrxns)):  # Change the '?' for the number of N-glycans
                if NGlyrxns[i] == 'XXX-M5-unfold-UBIQP_c + ? h2o_c + RAD23A_c =>  XXX-unfold-UBIQP-RAD23A_c + ? acgam_c + ? man_c':
                    h2o = str(7 * int(number_Nglycans))
                    acgam = str(2 * int(number_Nglycans))
                    man = str(5 * int(number_Nglycans))
                    NGlyrxns[i] = 'XXX-M5-unfold-UBIQP_c + ' + h2o + ' h2o_c + RAD23A_c =>  XXX-unfold-UBIQP_c-RAD23A_c + ' + acgam + ' acgam_c + ' + man + ' man_c'
                if '?' in NGlyrxns[i]:
                    if number_Nglycans == '1':
                        NGlyrxns[i] = NGlyrxns[i].replace('? ', '')
                    else:
                        NGlyrxns[i] = NGlyrxns[i].replace('?', number_Nglycans)

            for reaction in NGlyrxns:
                self.rxns.append(reaction)
            for reactionName in NGlyrxnNames:
                self.rxnNames.append(reactionName)

            self.connector = 'XXX-M8B_r'

        # Add COPII reactions
        if self.connector == 'XXX-M8B_r':
            copii_rxns = []
            copii_names = []
            copii_rxns, copii_names = self.addPathway('COPII_NG', copii_rxns, copii_names)
            for i in range(len(copii_rxns)):
                copii_rxns[i] = copii_rxns[i].replace("!", str(copii_coeff))
                self.rxns.append(copii_rxns[i])
                self.rxnNames.append(copii_names[i])

            self.connector = 'XXX-M3-GN2[g]'

            # Add if protein is EPO (Human)
            if self.entryID == 'P01588':
                self.rxns.append('P01588_c --> EPO_Human_c')
                self.rxnNames.append('make_EPO')
                self.rxns, self.rxnNames = self.addPathway('Golgi processing (EPO specific)', self.rxns, self.rxnNames)
                self.connector = 'XXX-M3-GN4-GL4-NA4-F[g]'
                self.entryID = 'EPO_Human'

        elif self.connector == 'XXX-dgpi_cho_r':
            copii_rxns = []
            copii_names = []
            copii_rxns, copii_names = self.addPathway('COPII_GPI', copii_rxns, copii_names)
            for i in range(len(copii_rxns)):
                copii_rxns[i] = copii_rxns[i].replace("!", str(copii_coeff))
                self.rxns.append(copii_rxns[i])
                self.rxnNames.append(copii_names[i])

            self.connector = 'XXX-dgpi_cho[g]'

        elif self.connector == 'XXX_DSB_r':
            self.rxns, self.rxnNames = self.addPathway('COPII_DSB', self.rxns, self.rxnNames, self.rxnComponents)

            copii_rxns = []
            copii_names = []
            copii_rxns, copii_names = self.addPathway('COPII-canonical', copii_rxns, copii_names, self.rxnComponents)
            for i in range(len(copii_rxns)):
                copii_rxns[i] = copii_rxns[i].replace("!", str(copii_coeff))
                self.rxns.append(copii_rxns[i])
                self.rxnNames.append(copii_names[i])

            self.connector = 'XXX[g]'

        elif self.connector == 'XXX_r':
            self.rxns, self.rxnNames = self.addPathway('COPII-normal', self.rxns, self.rxnNames)

            copii_rxns = []
            copii_names = []
            copii_rxns, copii_names = self.addPathway('COPII-canonical', copii_rxns, copii_names)
            for i in range(len(copii_rxns)):
                copii_rxns[i] = copii_rxns[i].replace("!", str(copii_coeff))
                self.rxns.append(copii_rxns[i])
                self.rxnNames.append(copii_names[i])

            self.connector = 'XXX[g]'

        # Add O-glycosylation reactions
        if self.OG != '0':  # Has O-glycans?
            self.rxns.append(self.connector + ' --> XXX_preOG_g')
            self.rxnNames.append('Start_OG')
            number_Oglycans = str(self.OG)  # Get number of O-Glycans
            OGlyrxns = []
            OGlyrxnNames = []

            OGlyrxns, OGlyrxnNames = self.addPathwayFromCondition('OG>0', OGlyrxns, OGlyrxnNames)
            for i in range(len(OGlyrxns)):  # Change the '?' for the number of O-glycans
                if '?' in OGlyrxns[i]:
                    if number_Oglycans == '1':
                        OGlyrxns[i] = OGlyrxns[i].replace('? ', '')
                    else:
                        OGlyrxns[i] = OGlyrxns[i].replace('?', number_Oglycans)

            for reaction in OGlyrxns:
                self.rxns.append(reaction)
            for reactionName in OGlyrxnNames:
                self.rxnNames.append(reactionName)

            self.connector = 'XXX-Core2_g'

        # Stay in Golgi if protein is localized there
        if self.subloc == '[g]' or self.subloc == '[gm]':
            self.location = self.subloc
            if 'SP_degradation' in self.rxnNames:
                SPaas = self.count_AAs(sequence[0:22])  # Amino acids in signal peptide assuming length is 22 on average
                self.rxns[self.rxnNames.index('SP_degradation')] = self.substitute_AAs_count(self.rxns[self.rxnNames.index('SP_degradation')], SPaas)
            if 'Ubiquitination_degradation' in self.rxnNames:
                new_aas = self.count_AAs(sequence[22:])
                self.rxns[self.rxnNames.index('Ubiquitination_degradation')] = self.rxns[self.rxnNames.index('Ubiquitination_degradation')].replace("!", "?")
                self.rxns[self.rxnNames.index('Ubiquitination_degradation')] = self.substitute_AAs_count(self.rxns[self.rxnNames.index('Ubiquitination_degradation')], new_aas)
                self.rxns[self.rxnNames.index('Ubiquitination_degradation')] = self.rxns[self.rxnNames.index('Ubiquitination_degradation')].replace("?", str(L-22))
            # Add GPRs
            self.GPRs = self.getGPRsFromRxnNames(self.rxnNames)
            # Change the reaction names and their formulas to include the UniProt Identifier
            for i in range(len(self.rxns)):
                self.rxns[i] = self.insert_prot_name_in_rxnFormula(self.rxns[i], self.entryID)
            for i in range(len(self.rxnNames)):
                self.rxnNames[i] = self.insert_prot_name_in_rxnName(self.rxnNames[i], self.entryID)
            self.rxns.append(self.entryID + self.location + ' --> ')
            self.rxnNames.append(self.entryID + '_Final_demand')
            self.GPRs.append('')
            self.rxns, self.rxnNames, self.GPRs = self.addCanonicalRxns(self.rxns, self.rxnNames, self.GPRs)
            return self.rxns, self.rxnNames, self.GPRs

        # Add COPI
        elif self.subloc == '[r]' or self.subloc == '[rm]':
            self.rxns.append(self.connector + ' --> XXX_preCOPI_g')
            self.rxnNames.append('Start_COPI')

            copi_rxns = []
            copi_names = []
            copi_rxns, copi_names = self.addPathway('COPI', copi_rxns, copi_names)
            for i in range(len(copi_rxns)):
                copi_rxns[i] = copi_rxns[i].replace("!", str(copi_coeff))
                self.rxns.append(copi_rxns[i])
                self.rxnNames.append(copi_names[i])

            self.connector = 'XXX_mature_r'
            self.location = self.subloc
            if self.location == '[r]':
                self.rxns.append(self.connector + ' --> ')
                self.rxnNames.append(self.entryID + '_Final_demand')
            elif self.location == '[rm]':
                self.rxns.append(self.connector + ' --> XXX_mature' + self.location)
                self.rxnNames.append('Final_location_' + self.location)
                self.rxns.append('XXX_mature' + self.location + ' --> ')
                self.rxnNames.append(self.entryID + '_Final_demand')

        # Add Clathrin vesicles
        elif self.subloc == '[x]' or self.subloc == '[l]' or self.subloc == '[d]':
            self.rxns.append(self.connector + ' --> XXX-preClathrin_g')
            self.rxnNames.append('Start_Clathrin_vesicle')

            clath_rxns = []
            clath_names = []
            clath_rxns, clath_names = self.addPathway('Clathrin vesicles', clath_rxns, clath_names)
            for i in range(len(clath_rxns)):
                clath_rxns[i] = clath_rxns[i].replace("!", str(self.clathrin_coeff))
                self.rxns.append(clath_rxns[i])
                self.rxnNames.append(clath_names[i])

            self.connector = 'XXX_mature_cv'
            self.location = self.subloc
            self.rxns.append(self.connector + ' --> XXX_mature' + self.location)
            self.rxnNames.append('Final_location_' + self.location)
            self.rxns.append('XXX_mature' + self.location + ' --> ')
            self.rxnNames.append(self.entryID + '_Final_demand')

        # Send to corresponding location
        else:
            self.location = self.subloc
            self.rxns.append(self.connector + ' --> XXX-preSV_g')
            self.rxnNames.append('Start_Secretion')

            sv_rxns = []
            sv_names = []
            sv_rxns, sv_names = self.addPathway('SV', sv_rxns, sv_names, self.rxnComponents)
            for i in range(len(sv_rxns)):
                sv_rxns[i] = sv_rxns[i].replace("!", str(self.clathrin_coeff))
                self.rxns.append(sv_rxns[i])
                self.rxnNames.append(sv_names[i])

            if self.location == '':
                self.location = '_e'
            self.rxns.append('XXX_mature_sv' + ' --> XXX_mature' + self.location)
            self.rxnNames.append('Final_location_' + self.location)
            self.rxns.append('XXX_mature' + self.location + ' --> ')
            self.rxnNames.append('Final_demand')

        # Add coefficients to SP_degradation and Ubiquitin_degradation reactions (if applicable)
        if 'SP_degradation' in self.rxnNames:
            SPaas = self.count_AAs(self.sequence[0:22])  # Amino acids in signal peptide assuming length is 22 on average
            self.rxns[self.rxnNames.index('SP_degradation')] = self.substitute_AAs_count(self.rxns[self.rxnNames.index('SP_degradation')], SPaas)
        if 'Ubiquitination_degradation' in self.rxnNames:
            new_aas = self.count_AAs(self.sequence[22:])
            self.rxns[self.rxnNames.index('Ubiquitination_degradation')] = self.rxns[self.rxnNames.index('Ubiquitination_degradation')].replace("!", "?")
            self.rxns[self.rxnNames.index('Ubiquitination_degradation')] = self.substitute_AAs_count(self.rxns[self.rxnNames.index('Ubiquitination_degradation')], new_aas)
            self.rxns[self.rxnNames.index('Ubiquitination_degradation')] = self.rxns[self.rxnNames.index('Ubiquitination_degradation')].replace("?", str(self.L-22))

        self.GPRs = self.getGPRsFromRxnNames(self.rxnNames)

        # Change the reaction names and their formulas to include the UniProt Identifier
        for i in range(len(self.rxns)):
            self.rxns[i] = self.insert_prot_name_in_rxnFormula(self.rxns[i], self.entryID)
        for i in range(len(self.rxnNames)):
            self.rxnNames[i] = self.insert_prot_name_in_rxnName(self.rxnNames[i], self.entryID)
        self.rxns, self.rxnNames, self.GPRs = self.addCanonicalRxns(self.rxns, self.rxnNames, self.GPRs)
        return self.rxns, self.rxnNames, self.GPRs

