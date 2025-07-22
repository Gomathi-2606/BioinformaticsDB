DROP DATABASE IF EXISTS bioinformaticsdb;
CREATE DATABASE IF NOT EXISTS BioinformaticsDB;
USE BioinformaticsDB;
-- Table-1
CREATE TABLE Researchers (
 researcher_id INT PRIMARY KEY AUTO_INCREMENT,
 name VARCHAR(100) NOT NULL,
 email VARCHAR(100) UNIQUE,
 phone_number VARCHAR(15),
 department VARCHAR(100),
 institution VARCHAR(150),
 country VARCHAR(100),
 date_joined DATE
);
INSERT INTO Researchers (name, email, phone_number, department, institution, country, date_joined) VALUES
('Dr. Alice Johnson', 'alice.johnson@example.com', '+1234567890', 'Genomics', 'BioTech Institute', 'USA', '2020-05-10'),
('Dr. Bob Smith', 'bob.smith@example.com', '+1987654321', 'Proteomics', 'LifeScience Labs', 'Canada', '2021-01-15'),
('Dr. Clara Lee', 'clara.lee@example.com', '+447700900123', 'Bioinformatics', 'Genome Center', 'UK', '2019-09-25'),
('Dr. David Kim', 'david.kim@example.com', '+8201051234567', 'Immunology', 'Seoul MedTech', 'South Korea', '2022-03-11'),
('Dr. Eva Müller', 'eva.mueller@example.com', '+4915123456789', 'Neuroscience', 'NeuroLab', 'Germany', '2018-07-30'),
('Dr. Farhan Ali', 'farhan.ali@example.com', '+911234567890', 'Microbiology', 'Indian Bio Inst.', 'India', '2020-12-05'),
('Dr. Grace Chen', 'grace.chen@example.com', '+8613800138000', 'Pharmacology', 'Shanghai Health Univ.', 'China', '2021-06-22'),
('Dr. Hiro Tanaka', 'hiro.tanaka@example.com', '+81312345678', 'Toxicology', 'Tokyo Research Inst.', 'Japan', '2017-11-03'),
('Dr. Isabella Rossi', 'isabella.rossi@example.com', '+390612345678', 'Pathology', 'MedScience Rome', 'Italy', '2019-04-18'),
('Dr. Juan Pérez', 'juan.perez@example.com', '+5215512345678', 'Virology', 'Mexican National Lab', 'Mexico', '2023-02-10');

-- Table-2
CREATE TABLE Genes (
 gene_id INT PRIMARY KEY AUTO_INCREMENT,
 gene_name VARCHAR(100) NOT NULL,
 description TEXT,
 organism VARCHAR(100),
 chromosome_location VARCHAR(100),
 start_position BIGINT,
 end_position BIGINT,
 strand ENUM('+', '-'),
 gene_type VARCHAR(50),
 gene_sequence TEXT,
 discovery_date DATE
);
INSERT INTO Genes (gene_name, description, organism, chromosome_location,start_position, end_position, strand, gene_type, gene_sequence, discovery_date) VALUES
('BRCA1', 'Breast cancer type 1 susceptibility protein', 'Homo sapiens', '17q21.31',43044295, 43170245, '+', 'protein_coding', 'ATGCGTACGTTAGC...', '1994-10-01'),
('TP53', 'Tumor protein p53', 'Homo sapiens', '17p13.1',7668402, 7687550, '-', 'protein_coding', 'ATGCCAGTCCTAGT...', '1979-05-20'),
('EGFR', 'Epidermal growth factor receptor', 'Homo sapiens', '7p11.2',55086714, 55279321, '+', 'protein_coding', 'ATGGAGCTGCCCGT...', '1984-11-12'),
('MYC', 'MYC proto-oncogene', 'Homo sapiens', '8q24.21',128748315, 128753680, '+', 'protein_coding', 'ATGGATTTTTTCCG...', '1983-06-15'),
('CFTR', 'Cystic fibrosis transmembrane conductance regulator', 'Homo sapiens', '7q31.2',117120016, 117308718, '-', 'protein_coding', 'ATGAGAGCGAGACT...', '1989-09-08'),
('APOE', 'Apolipoprotein E', 'Homo sapiens', '19q13.32',44905754, 44909393, '-', 'protein_coding', 'ATGCCGATGCTGGA...', '1985-03-01'),
('HBB', 'Hemoglobin subunit beta', 'Homo sapiens', '11p15.4',5225466, 5227901, '+', 'protein_coding', 'ATGGTGCACCTGAC...', '1960-12-15'),
('FMR1', 'Fragile X mental retardation 1', 'Homo sapiens', 'Xq27.3',146993207, 147032564, '+', 'protein_coding', 'ATGGAGGCCGGGAC...', '1991-03-29'),
('DMD', 'Dystrophin', 'Homo sapiens', 'Xp21.2',3110000, 3333942, '-', 'protein_coding', 'ATGGCAGCTGCTGA...', '1987-07-10'),
('HTT', 'Huntingtin', 'Homo sapiens', '4p16.3',3074877, 3245490, '+', 'protein_coding', 'ATGGCGACCCTGGG...', '1993-03-24');

-- Table-3
CREATE TABLE Proteins (
    protein_id INT PRIMARY KEY AUTO_INCREMENT,
    gene_id INT,
    protein_name VARCHAR(100),
    sequence TEXT,
    amino_acid_length INT,
    molecular_weight DECIMAL(10,3),
    isoelectric_point DECIMAL(5,2),
    cellular_location VARCHAR(100),
    protein_function TEXT,
    structure_3d_url VARCHAR(255),
    FOREIGN KEY (gene_id) REFERENCES Genes(gene_id) ON DELETE CASCADE
);
INSERT INTO Proteins (gene_id, protein_name, sequence, amino_acid_length,molecular_weight, isoelectric_point, cellular_location, protein_function, structure_3d_url) VALUES
(1, 'Breast cancer type 1 susceptibility protein', 'MSEQDTVQKL...', 1863, 220.000, 5.80, 'Nucleus','DNA repair', 'https://proteinstructures.org/brca1'),
(2, 'Cellular tumor antigen p53', 'MEEPQSDPSV...', 393, 43.700, 6.50, 'Nucleus','Tumor suppression', 'https://proteinstructures.org/tp53'),
(3, 'Epidermal growth factor receptor', 'MRPSGTAGAALL...', 1210, 134.000, 6.20, 'Cell membrane','Cell signaling', 'https://proteinstructures.org/egfr'),
(4, 'Myc proto-oncogene protein', 'MPLNVSFTNRN...', 439, 48.000, 5.50, 'Nucleus','Transcription regulation', 'https://proteinstructures.org/myc'),
(5, 'Cystic fibrosis transmembrane conductance regulator', 'MSDQEQPVQVTF...', 1480, 168.000, 6.80, 'Plasma membrane','Chloride ion transport', 'https://proteinstructures.org/cftr'),
(6, 'Apolipoprotein E', 'MKAFLCLWLLL...', 317, 36.000, 5.90, 'Extracellular','Lipid transport', 'https://proteinstructures.org/apoe'),
(7, 'Hemoglobin subunit beta', 'MVHLTPEEKSAV...', 147, 15.900, 6.80, 'Red blood cell','Oxygen transport', 'https://proteinstructures.org/hbb'),
(8, 'Fragile X mental retardation protein 1', 'MDQLSPQLVQ...', 632, 70.000, 8.20, 'Cytoplasm','RNA binding', 'https://proteinstructures.org/fmr1'),
(9, 'Dystrophin', 'MSKAEIEQEAL...', 3685, 427.000, 5.10, 'Cytoskeleton','Muscle fiber stability', 'https://proteinstructures.org/dmd'),
(10, 'Huntingtin', 'MENLVQFTPAF...', 3144, 350.000, 5.30, 'Cytoplasm','Neuronal signaling', 'https://proteinstructures.org/htt');

-- Table-4
CREATE TABLE Samples (
 sample_id INT PRIMARY KEY AUTO_INCREMENT,
 sample_name VARCHAR(100),
 sample_type VARCHAR(50),
 tissue_type VARCHAR(100),
 collected_by VARCHAR(100),
 researcher_id INT,
 collection_date DATE,
 organism VARCHAR(100),
 storage_location VARCHAR(100),
 storage_temperature DECIMAL(4, 1),
 quality_score INT,
 FOREIGN KEY (researcher_id) REFERENCES Researchers(researcher_id) ON DELETE SET
NULL
);
INSERT INTO Samples (sample_name, sample_type, tissue_type, collected_by, researcher_id,collection_date, organism, storage_location, storage_temperature, quality_score) VALUES
('Sample A1', 'Tissue', 'Breast tissue', 'Dr. Alice Johnson', 1,'2021-01-15', 'Homo sapiens', 'Freezer A - Rack 2', -80.0, 95),
('Sample B2', 'Blood', 'Whole blood', 'Dr. Bob Smith', 2,'2021-02-10', 'Homo sapiens', 'Freezer B - Rack 1', -70.0, 92),
('Sample C3', 'Saliva', 'Buccal swab', 'Dr. Clara Lee', 3,'2020-12-05', 'Homo sapiens', 'Fridge 1 - Shelf 3', 4.0, 90),
('Sample D4', 'Tissue', 'Liver tissue', 'Dr. David Kim', 4,'2022-03-22', 'Homo sapiens', 'Freezer C - Rack 5', -80.0, 93),
('Sample E5', 'Plasma', 'Blood plasma', 'Dr. Eva Müller', 5,'2019-11-08', 'Homo sapiens', 'Freezer A - Rack 3', -75.0, 91),
('Sample F6', 'Tissue', 'Colon tissue', 'Dr. Farhan Ali', 6,'2021-09-18', 'Homo sapiens', 'Freezer D - Rack 1', -80.0, 94),
('Sample G7', 'Urine', 'Urine', 'Dr. Grace Chen', 7,'2020-06-25', 'Homo sapiens', 'Fridge 2 - Shelf 1', 4.0, 88),
('Sample H8', 'CSF', 'Cerebrospinal fluid', 'Dr. Hiro Tanaka', 8,'2023-01-14', 'Homo sapiens', 'Freezer E - Rack 4', -78.0, 89),
('Sample I9', 'Tissue', 'Brain tissue', 'Dr. Isabella Rossi', 9,'2018-08-30', 'Homo sapiens', 'Freezer B - Rack 6', -82.0, 96),
('Sample J10', 'Serum', 'Blood serum', 'Dr. Juan Pérez', 10,'2023-05-11', 'Homo sapiens', 'Freezer F - Rack 2', -70.0, 93);

-- Table-5
CREATE TABLE Experiments (
 experiment_id INT PRIMARY KEY AUTO_INCREMENT,
 sample_id INT,
 researcher_id INT,
 experiment_name VARCHAR(100),
 experiment_type VARCHAR(50),
 technique_used VARCHAR(100),
 date DATE,
 status ENUM('Pending', 'Completed', 'Failed') DEFAULT 'Pending',
 results TEXT,
 notes TEXT,
 FOREIGN KEY (sample_id) REFERENCES Samples(sample_id) ON DELETE CASCADE,
 FOREIGN KEY (researcher_id) REFERENCES Researchers(researcher_id) ON DELETE SET
NULL
);
INSERT INTO Experiments (sample_id, researcher_id, experiment_name, experiment_type, technique_used,date, status, results, notes) VALUES
(1, 1, 'Gene Expression Analysis', 'RNA-Seq', 'Illumina HiSeq','2021-02-20', 'Completed', 'Expression levels quantified successfully.', 'No anomalies detected.'),
(2, 2, 'Protein Profiling', 'Proteomics', 'LC-MS/MS','2021-03-15', 'Completed', 'Protein expression profile generated.', 'Data is clean and complete.'),
(3, 3, 'Methylation Study', 'Epigenetics', 'Bisulfite Sequencing','2021-01-10', 'Completed', 'Methylation patterns identified.', 'Results match expected tissue patterns.'),
(4, 4, 'Transcriptome Assembly', 'RNA-Seq', 'PacBio IsoSeq','2022-04-05', 'In Progress', 'Partial transcriptome assembled.', 'Awaiting further read coverage.'),
(5, 5, 'SNP Genotyping', 'Genotyping', 'Affymetrix Array','2020-12-22', 'Completed', 'Several SNPs identified.', 'SNP distribution normal.'),
(6, 6, 'Gene Knockout Validation', 'CRISPR', 'qPCR','2021-10-03', 'Completed', 'Gene knockout confirmed by qPCR.', 'Target gene expression absent.'),
(7, 7, 'Metabolite Analysis', 'Metabolomics', 'GC-MS','2020-07-29', 'Completed', 'Metabolite profiles generated.', 'Sample preparation was optimal.'),
(8, 8, 'Protein Structure Prediction', 'Structural Biology', 'AlphaFold','2023-02-18', 'Completed', '3D structure model predicted successfully.', 'Model meets confidence thresholds.'),
(9, 9, 'Neurotransmitter Quantification', 'Neurobiology', 'HPLC','2019-09-05', 'Completed', 'Neurotransmitters quantified.', 'All results within expected ranges.'),
(10, 10, 'Immunoassay Screening', 'Immunology', 'ELISA','2023-06-10', 'In Progress', 'Preliminary data obtained.', 'Requires duplicate verification.');

-- Table-6
CREATE TABLE GeneAnnotations (
 annotation_id INT PRIMARY KEY AUTO_INCREMENT,
 gene_id INT,
 annotation TEXT,
 ontology_term VARCHAR(100),
 source VARCHAR(100),
 confidence_score DECIMAL(4,2),
 date_added DATE,
 added_by VARCHAR(100),
 FOREIGN KEY (gene_id) REFERENCES Genes(gene_id) ON DELETE CASCADE
);
INSERT INTO GeneAnnotations (gene_id, annotation, ontology_term, source, confidence_score,date_added, added_by) VALUES
(1, 'Involved in DNA repair pathways.', 'GO:0006281', 'UniProt', 0.95, '2021-03-01', 'Dr. Alice Johnson'),
(2, 'Regulates cell cycle.', 'GO:0007049', 'Ensembl', 0.87, '2021-04-15', 'Dr. Bob Lee'),
(3, 'Plays role in immune response.', 'GO:0006955', 'NCBI', 0.90, '2021-06-10', 'Dr. Clara Zhang'),
(4, 'Acts in homologous recombination repair.', 'GO:0035825', 'UniProt', 0.91, '2021-07-20', 'Dr. Alice Johnson'),
(5, 'Cell cycle checkpoint control.', 'GO:0072395', 'GeneCards', 0.85, '2021-08-05', 'Dr. Bob Lee'),
(6, 'Oncogene with transcription regulation activity.', 'GO:0003700', 'UniProt', 0.88, '2021-09-12', 'Dr. David Kim'),
(7, 'Chloride ion transport protein.', 'GO:0005254', 'NCBI', 0.92, '2020-10-30', 'Dr. Eva Müller'),
(8, 'Involved in lipid metabolism.', 'GO:0006629', 'GeneCards', 0.89, '2021-01-20', 'Dr. Farhan Ali'),
(9, 'Hemoglobin beta chain binds oxygen.', 'GO:0015671', 'UniProt', 0.96, '2021-02-18', 'Dr. Grace Chen'),
(10, 'RNA-binding involved in neuronal development.', 'GO:0003723', 'Ensembl', 0.84, '2022-01-05', 'Dr. Hiro Tanaka');

-- Table-7
CREATE TABLE Diseases (
 disease_id INT PRIMARY KEY AUTO_INCREMENT,
 disease_name VARCHAR(100),
 gene_id INT,
 disease_type VARCHAR(100),
 affected_population VARCHAR(100),
 inheritance_pattern VARCHAR(100),
 prevalence_rate DECIMAL(5,2),
 description TEXT,
 FOREIGN KEY (gene_id) REFERENCES Genes(gene_id) ON DELETE SET NULL
);
INSERT INTO Diseases (disease_name, gene_id, disease_type, affected_population,inheritance_pattern, prevalence_rate, description) VALUES
('Hereditary Breast and Ovarian Cancer', 1, 'Cancer','Females of European descent', 'Autosomal Dominant', 0.20,'Inherited mutations in BRCA1 significantly increase cancer risk.'),
('Li-Fraumeni Syndrome', 2, 'Cancer','Global', 'Autosomal Dominant', 0.01,'Caused by mutations in TP53 gene, leading to a variety of cancers.'),
('Lung Adenocarcinoma', 3, 'Cancer','Asian non-smokers', 'Somatic Mutation', 0.18,'EGFR mutations are commonly observed in non-smoking lung cancer patients.'),
('Burkitt Lymphoma', 4, 'Cancer','African children', 'Translocation', 0.05,'MYC translocation is associated with aggressive B-cell lymphoma.'),
('Cystic Fibrosis', 5, 'Genetic Disorder','Caucasian population', 'Autosomal Recessive', 0.04,'Mutations in CFTR gene lead to thick mucus and respiratory issues.'),
('Alzheimer’s Disease (Late Onset)', 6, 'Neurodegenerative','Elderly population', 'Polygenic', 0.12,'APOE ε4 allele increases risk for late-onset Alzheimer’s.'),
('Beta Thalassemia', 7, 'Blood Disorder','Mediterranean and South Asian populations', 'Autosomal Recessive', 0.10,'HBB gene mutations reduce or eliminate production of beta globin.'),
('Fragile X Syndrome', 8, 'Genetic Disorder','Males globally', 'X-linked Dominant', 0.02,'Expansion of CGG repeats in FMR1 gene causes intellectual disability.'),
('Duchenne Muscular Dystrophy', 9, 'Muscular Disorder','Male children', 'X-linked Recessive', 0.015,'Caused by mutations in the DMD gene affecting dystrophin production.'),
('Huntington’s Disease', 10, 'Neurodegenerative','European ancestry', 'Autosomal Dominant', 0.005,'HTT gene CAG repeat expansion causes progressive brain disorder.');

-- Table-8
CREATE TABLE Pathways (
 pathway_id INT PRIMARY KEY AUTO_INCREMENT,
 pathway_name VARCHAR(100),
 pathway_type VARCHAR(100),
 description TEXT,
 source_database VARCHAR(100),
 created_by VARCHAR(100),
 date_created DATE
);
INSERT INTO Pathways (pathway_name, pathway_type, description, source_database,created_by, date_created) VALUES
('Homologous Recombination', 'DNA Repair','Repair of double-strand breaks via homologous recombination.','KEGG', 'Dr. Alice Johnson', '2021-04-15'),
('Cell Cycle Control', 'Cell Cycle','Regulates progression through phases of the cell cycle.','Reactome', 'Dr. Bob Lee', '2021-05-10'),
('EGFR Signaling Pathway', 'Signal Transduction','Pathway activated by epidermal growth factor receptor binding.','KEGG', 'Dr. Clara Lee', '2021-06-02'),
('MYC Activation Pathway', 'Transcription Regulation','Controls expression of genes promoting cell growth and division.','WikiPathways', 'Dr. David Kim', '2021-07-01'),
('Cystic Fibrosis Transmembrane Conductance', 'Ion Transport','Regulates chloride ion flow across membranes in epithelial cells.','KEGG', 'Dr. Eva Müller', '2020-09-15'),
('Apolipoprotein E Pathway', 'Lipid Metabolism','Involved in cholesterol transport and uptake in the brain.','Reactome', 'Dr. Farhan Ali', '2021-03-18'),
('Hemoglobin Function Pathway', 'Oxygen Transport','Facilitates oxygen transport from lungs to tissues.','NCBI BioSystems', 'Dr. Grace Chen', '2021-01-22'),
('Fragile X Pathway', 'RNA Processing','Pathway affected by loss of FMRP, involved in RNA binding.','KEGG', 'Dr. Hiro Tanaka', '2022-01-04'),
('Dystrophin-Glycoprotein Complex', 'Muscle Integrity','Stabilizes muscle fiber membrane during contraction.','Reactome', 'Dr. Isabella Rossi', '2020-12-11'),
('Huntingtin Pathway', 'Neurodegeneration','Involves disrupted cellular processes in Huntington\'s disease.','WikiPathways', 'Dr. Juan Pérez', '2023-03-01');

-- Table-9
CREATE TABLE GenePathways (
 gene_id INT,
 pathway_id INT,
 role_in_pathway VARCHAR(100),
 PRIMARY KEY (gene_id, pathway_id),
 FOREIGN KEY (gene_id) REFERENCES Genes(gene_id) ON DELETE CASCADE,
 FOREIGN KEY (pathway_id) REFERENCES Pathways(pathway_id) ON DELETE CASCADE
);
INSERT INTO GenePathways (gene_id, pathway_id, role_in_pathway) VALUES
(1, 1, 'Core component of repair machinery'),
(2, 2, 'Regulates checkpoint transitions'),
(3, 3, 'Activates downstream signaling cascade'),
(4, 4, 'Drives transcription of growth-related genes'),
(5, 5, 'Chloride ion channel protein'),
(6, 6, 'Mediates lipid transport and receptor binding'),
(7, 7, 'Binds and transports oxygen in red blood cells'),
(8, 8, 'RNA-binding involved in synaptic function'),
(9, 9, 'Links cytoskeleton to extracellular matrix'),
(10, 10, 'Modulates neuronal survival and transcription');

-- Table-10
CREATE TABLE Publications (
 publication_id INT PRIMARY KEY AUTO_INCREMENT,
 title VARCHAR(200),
 journal VARCHAR(100),
 authors TEXT,
 doi VARCHAR(100),
 year INT,
 gene_id INT,
 researcher_id INT,
 abstract TEXT,
 publication_date DATE,
 FOREIGN KEY (gene_id) REFERENCES Genes(gene_id) ON DELETE SET NULL,
 FOREIGN KEY (researcher_id) REFERENCES Researchers(researcher_id) ON DELETE SET
NULL
);
INSERT INTO Publications (title, journal, authors, doi, year,gene_id, researcher_id, abstract, publication_date) VALUES
('Role of BRCA1 in DNA Repair', 'Nature Genetics', 'Alice Johnson, Bob Smith','10.1038/ng1234', 2021, 1, 1,'This study elucidates the mechanism of BRCA1-mediated homologous recombination.','2021-06-10'),
('TP53 as a Master Regulator of Apoptosis', 'Cell', 'Bob Lee, Clara Zhang','10.1016/j.cell.2021.02.015', 2021, 2, 2,'A comprehensive analysis of TP53’s role in cell cycle arrest and apoptosis.','2021-05-18'),
('EGFR Mutations in Lung Adenocarcinoma', 'Journal of Clinical Oncology', 'Clara Lee, David Kim','10.1200/JCO.21.00321', 2021, 3, 3,'Findings highlight the prevalence and treatment implications of EGFR mutations.','2021-07-01'),
('MYC-Driven Tumorigenesis Pathways', 'Cancer Cell', 'David Kim, Eva Müller','10.1016/j.ccell.2021.04.012', 2021, 4, 4,'Study explores transcriptional targets of MYC in cancer proliferation.','2021-08-03'),
('CFTR Gene Mutations and Cystic Fibrosis', 'Nature Medicine', 'Eva Müller, Farhan Ali','10.1038/nm12345', 2020, 5, 5,'A detailed genetic and clinical study of CFTR mutations and their phenotypes.','2020-12-20'),
('APOE ε4 Allele in Alzheimer’s Disease', 'Neuron', 'Farhan Ali, Grace Chen','10.1016/j.neuron.2021.06.004', 2021, 6, 6,'The APOE ε4 allele is correlated with early-onset Alzheimer’s risk.','2021-11-02'),
('Beta Globin Variants in Thalassemia', 'Blood', 'Grace Chen, Hiro Tanaka','10.1182/blood.2021.88.5.1234', 2021, 7, 7,'Beta globin gene mutations and their effect on thalassemia severity.','2021-09-14'),
('FMR1 Gene and Fragile X Syndrome Mechanisms', 'Molecular Psychiatry', 'Hiro Tanaka, Isabella Rossi','10.1038/mp.2021.76', 2022, 8, 8,'Insights into FMRP functions and CGG repeat expansions.','2022-01-22'),
('Dystrophin Gene Deletions in DMD', 'Human Molecular Genetics', 'Isabella Rossi, Juan Pérez','10.1093/hmg/ddab321', 2020, 9, 9,'Detailed mapping of dystrophin deletions and clinical correlations.','2020-10-30'),
('Huntingtin Protein in Neurodegeneration', 'Brain', 'Juan Pérez, Alice Johnson','10.1093/brain/awab101', 2023, 10, 10,'HTT gene’s role in neurodegenerative processes and potential therapeutic targets.','2023-03-05');

-- Table-11
CREATE TABLE Mutations (
 mutation_id INT PRIMARY KEY AUTO_INCREMENT,
 gene_id INT,
 mutation_type VARCHAR(100),
 nucleotide_change VARCHAR(50),
 protein_change VARCHAR(50),
 clinical_significance VARCHAR(100),
 detection_method VARCHAR(100),
 date_reported DATE,
 reference_study VARCHAR(255),
 FOREIGN KEY (gene_id) REFERENCES Genes(gene_id) ON DELETE CASCADE
);
INSERT INTO Mutations (gene_id, mutation_type, nucleotide_change, protein_change,clinical_significance, detection_method, date_reported, reference_study) VALUES
(1, 'Frameshift', 'c.5266dupC', 'p.Gln1756Profs', 'Pathogenic', 'Sanger Sequencing','2021-07-05', 'Nature Genetics, 2021'),
(2, 'Nonsense', 'c.637C>T', 'p.Arg213Ter', 'Pathogenic', 'NGS Panel','2021-05-10', 'Cell, 2021'),
(3, 'Missense', 'c.2573T>G', 'p.Leu858Arg', 'Likely Pathogenic', 'Whole Exome Sequencing','2021-06-15', 'J Clin Oncol, 2021'),
(4, 'Amplification', 'n/a', 'n/a', 'Oncogenic', 'FISH','2021-09-01', 'Cancer Cell, 2021'),
(5, 'Splice Site', 'c.2988+1G>A', 'p.Thr1030_Ser1031ins', 'Pathogenic', 'RT-PCR','2020-11-25', 'Nature Medicine, 2020'),
(6, 'Missense', 'c.388T>C', 'p.Cys130Arg', 'Risk Factor', 'SNP Genotyping','2021-02-08', 'Neuron, 2021'),
(7, 'Deletion', 'c.92+5_92+8del', 'p.Val31Glufs', 'Pathogenic', 'MLPA','2021-08-13', 'Blood, 2021'),
(8, 'Expansion', 'CGG>200 repeats', 'n/a', 'Pathogenic', 'Southern Blot','2022-01-20', 'Mol Psychiatry, 2022'),
(9, 'Large Deletion', 'Exons 45–50 del', 'p.Trp1469_Glu1676del', 'Pathogenic', 'MLPA','2020-10-05', 'Hum Mol Genet, 2020'),
(10, 'Trinucleotide Repeat', 'c.52CAG>n', 'p.Gln18poly', 'Pathogenic', 'Triplet-Primed PCR','2023-03-10', 'Brain, 2023');

SELECT*FROM Researchers;
SELECT organism,chromosome_location FROM Genes;
SELECT*FROM Mutations WHERE date_reported >2020;
SELECT*FROM Publications WHERE authors Like 'A%';
SELECT pathway_type, COUNT(*) AS total_pathways
FROM Pathways
GROUP BY pathway_type
HAVING COUNT(*) > 1;
SELECT 
    gp.gene_id,
    gp.pathway_id,
    gp.role_in_pathway,
    p.pathway_name,
    p.pathway_type
FROM 
    GenePathways gp
INNER JOIN 
    Pathways p ON gp.pathway_id = p.pathway_id;
SELECT 
    d.disease_name,
    d.disease_type,
    d.affected_population,
    d.inheritance_pattern,
    d.prevalence_rate,
    g.gene_name,
    g.chromosome_location
FROM 
    Diseases d
LEFT JOIN 
    Genes g ON d.gene_id = g.gene_id;
SELECT 
    g.gene_id,
    g.gene_name,
    g.chromosome_location,
    ga.annotation,
    ga.ontology_term,
    ga.source,
    ga.confidence_score
FROM 
    GeneAnnotations ga
RIGHT JOIN 
    Genes g ON ga.gene_id = g.gene_id;
SELECT COUNT(*) AS total_experiments
FROM Experiments;
SELECT AVG(quality_score) AS average_quality_score
FROM Samples;
SELECT 
    MAX(molecular_weight) AS max_molecular_weight, 
    MIN(molecular_weight) AS min_molecular_weight,
    MAX(amino_acid_length) AS max_amino_acid_length, 
    MIN(amino_acid_length) AS min_amino_acid_length,
    MAX(isoelectric_point) AS max_isoelectric_point, 
    MIN(isoelectric_point) AS min_isoelectric_point
FROM Proteins;












