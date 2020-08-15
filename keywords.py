from flashtext import KeywordProcessor
import re

# General logic of keywords in lists and dictionaries below:
	# To match people on as many generic terms as possible that encompass vast scientific areas and fields, 
	# but then to get bonus points for super specific overlaps in terminology

# Note: key words all have equal weight. 


# Keyword Lists & Dictionaries

# Singular Keywords (both super vague and super specific)
science = ["applied", "chemical", "theoretical", "electrical", "computer", "mechanical", "materials", "biomedical",
	"structural", "clinical", "global warming", "radiation", "observational", "translational",
	"computational", "open science", "inorganic", "surgical", "behavioral", "vascular", "quantiative",
	"airborne", "marijuana", "ketogenic", "universe", "engineering", "atmospheric", "stress",
	"biometrics", "EEG", "EKG", "myocardial ischemia-reperfusion", "skeletal", "operations",
	"fluorescence-based", "mitochondrial", "multiagent", "poaching", "radiology",
	"synthetic", "condensed matter", "organic", "metallurgy", "proteomics", "e. coli",
	"tropical", "hippocampus", "biomaterial analysis", "gravitational waves", "commutative",
	"high performance computing", "electronic", "open source", "Riemannian Hypothesis",
	"materials science", "synthesis", "resistance", "fungal", "bedside", "energy"
	"pulmonary inflammation", "biophysics", "networks", "structural", "physical",
	"biodegradable polymer", "biodegradable", "Csx3", "x-ray crystallography",  
	"HEP", "renewable energy", "solar", "radio", "pure", "yeast", "mammalian", 
	"Renal Fibrosis", "postnatal murine cardiac fibroblasts", "bacteriophage lysis", 
	"non-arteritic anterior ischemic optic neuropathy", "NAION", "biomechanics", "nonlinear",
	"lumbar", "aging", "degeneration", "alternating magnetic fields", "hydraulics", 
	"vegetated flow patches", "eutrophication", "anthropogenic", "photonics", "particle",  
	"linguistics", "adaptive", "oxygenation", "blood", "mechanics", "remote sensing"
	"nonlinearly coupled masses", "hall effect", "control", "differential equations",
	"mental", "feminism", "diversity", "endothelium", "migratory signaling", "p53", "S100B", 
	"physicoceanal ography", "climate change", "archaic myosin", "concussive",
	"Hypertrophy Cardiomyopathy", "abstract", "representation theory", "bromodomain", 
	"viral hemorrhagic fevers", "ketamine", "non-equilibrium phases", "AMO", "histone modifications", 
	"asymmetric synthesis", "neuroprotectant", "organic synthesis", "epigenetic aberrations",
	"Mukaiyama aldol addition reactions", "industrial", "KSHV", "optics", "diamond", "fiber",
	"silicon carbide", "sensing", "cell structure", "morphometrics", "learning theory", "aggregation",
	"spatial", "herpetofauna", "mismatch repair", "neutrinos", "stochastic control techniques",
	"phytoplankton", "tissue", "shrooms", "metabolic engineering", "modeling", "methodology", "mechanism" 
	"bio-interfaces", "climate", "methane flux", "bioacoustics", "clam", "oyster", "drug", "repetitive", 
	"ecophysiology", "climatological niche modeling", "Rheumatoid arthritis", "Graft Versus Host Disease", 
	"reproductive", "modulate", "preterm birth", "nuclear", "impulsivity", "primary cilia",
	"thioredoxin", "oxidative stress", "computational science", "biomodulatory materials", "bone tissue", 
	"regenerative", "MRSA", "coastal", "particle", "neutrino", "zoonotic", "release kinetics", 
	"chromatin", "antisense oligonucleotide", "dots", "regeneration", "deformity",
	"protein dynamics", "noisy intermediate-scale", "NISQ", "civil" 
	"hematopoiesis", "high energy", "thermodynamics", "biomaterials", "peripheral nerve", "brain", 
	"amygdala", "olfactory", "sensory", "prefrontal cortex", "somatic", "memory", "learning",  
	"reactive oxygen species", "damage repair", "physical", "medical sciences", "pre-med",
	"allergic inflammation", "patient-data-based", "ultracold", "petrology", "experimental", 
	"froth flotation", "antimicrobial", "meteorology", "bioremediation", "research scientist",
	"medicolegal", "endangered species", "pain", "diagnostic tools", "detector", "accelerator",
	"psychopathological", "images", "imaging", "radioactive", "topological materials",
	"forest", "alpine", "pond", "optical", "metalloenzymatic", "spectroscopy", "spatial analysis"
	"hardware", "software", "mass", "spectrometry", "synaptic reorganization", "liver", "microscopy",
	"infrastructure", "enteric", "nervous system", "biophysical", "apnea", "sleep", "energy storage",
	"orthopaedic", "algae", "water", "illness", "organometallic", "membrane technology", "metastasis",
	"suprachiasmatic nucleus", "hemmorhagic", "depressive", "Langerhans cell histiocytosis", 
	"fluid mechanics", "dynamics", "physiology", "ergonomics", "plasma", "norepinephrine",
	"fusion", "psilocybin", "black holes", "nitrogen", "oxygen", "computer vision", "nuclear imaging",
	"categorization", "craniofacial", "topology", "skin", "thyroid", "vector-borne", 
	"precision", "onshore", "offshore", "wind turbine", "subterranean", "battery", "folding", 
	"computational plasma", "aerogels", "adaptive radiations", "sensor deployment", "string theory", 
	"formation", "extragalactic", "functional", "nucleation", "aerosol", "fuels", "manganese", 
	"developmental", "bioengineering", "friction stir welding", "nanoindentation", "refractory",
	"non-human", "aquatic", "disinfection", "gold", "peptoid-functionalized", "mycorrhizal",
	"electronics", "pollution", "air", "dosimeter", "footprinting", "limbic system", 
	"biodiversity", "catalysis", "macromolecular", "supramolecular", "embedded"]

# Allows for variations of terms, associations between similar terms, 
# and marking terms under the same field umbrella through tuple keys
science_dict = { 
	# Neuroscience
	"neuro*": ["neuroscience", "neurology", "neurological", "neural"],
	("neuro*", "cognitive"): ["cognitive"],
	("neuro*", "neurodevelopment"): ["neurodevelopment", "neurodevelopmental"],
	("neuro*", "neuropharmacology"): ["neuropharmacology"],
	("neuro*", "neuroengineering"): ["neuroengineering", "neuro-engineering"],
	("neuro*", "neuropsychiatry"): ["neuropsychiatry"],
	("neuro*", "Alzheimer’s"): ["Alzheimer’s"],
	("neuro*", "ALS"): ["ALS", "Lou Gehrig's"],
	("neuro*", "Parkinson's"): ["Parkinson's"],
	("neuro*", "Huntington's"): ["Huntington's"],
	("neuro*", "Multiple Sclerosis"): ["Multiple Sclerosis"],
	("neuro*", "neurodegenerative"): ["neurodegenerative"],
	("neuro*", "neurobiology"): ["neurobiology"],
	("neuro*", "neuroinformatics"): ["neuroinformatics"],
	("neuro*", "neurosurgeon"): ["neurosurgeon"],
	("neuro*", "blood-brain barrier"): ["blood-brain barrier"],
	("neuro*", "neuropsychiatric"): ["neuropsychiatric"],
	("neuro*", "neurotrauma"): ["neurotrauma"],
	("neuro*", "neuroendocrinology"): ["neuroendocrine", "neuroendocrinology"],

	# Hormones
	"hormones*": ["hormones", "hormone"],
	("hormones*", "estrogen"): ["estrogen"],
	("hormones*", "endocrinology"): ["endocrinology", "endocrinologist"],
	("hormones*", "testosterone"): ["testosterone"],
	("hormones*", "endocrine"): ["endocrine"],
	("hormones*", "neuroendocrinology"): ["neuroendocrine", "neuroendocrinology"],

	# Diseases
	"disease*": ["disease", "diseases"],
	("disease*", "pathogens"): ["pathogens", "pathogen", "pathogenesis", "host-pathogen", "pathogenic"],
	("disease*", "pathology"): ["pathology", "pathologist"],
	("disease*", "paleopathology"): ["paleopathology", "paleopathologist"],
	("disease*", "virus"): ["virus", "viruses", "viral"],
	("disease*", "virology"): ["virology", "virologist"],
	("disease*", "ZIKA"): ["ZIKA"],
	("disease*", "epidemiology"): ["epidemiology", "epidemiological"],
	("disease*", "dengue"): ["dengue"],
	("disease*", "HIV/AIDS"): ["HIV", "AIDS"],
	("disease*", "coronavirus"): ["coronavirus", "covid-19", "coronviruses", "MERS-CoV", "MERS", "SARS-CoV", "SARS"],
	("disease*", "diabetes"): ["diabetes", "diabetic", "diabetics"],
	("disease*", "neurodegenerative"): ["neurodegenerative"],
	("disease*", "flu"): ["flu", "influenza"],
	("disease*", "Mad Cow Disease"): ["Mad Cow Disease"],
	("disease*", "Sindbis"): ["Sindbis"],
	("disease*", "Chikungunya"): ["Chikungunya "],

	# Disorders
	"disorder*": ["disorders", "disorder"],
	("disorder*", "autism"): ["autism", "autistic"],
	("disorder*", "epilepsy"): ["epilepsy", "epileptic"],
	("disorder*", "OCD"): ["OCD", "obsessive compulsive"],
	("disorder*", "schizophrenia"): ["schizophrenia", "schizophrenic"],
	("disorder*", "ADHD/ADD"): ["ADHD", "attention deficit hyperactivity", "ADD", "attention deficit"],
	("disorder*", "addiction"): ["addiction"],
	("disorder*", "movement disorders"): ["movement disorders", "movement disorder"],

	# Biology
	"biology*": ["biology", "biological", "biologist"], 
	("biology*", "enzymes*"): ["enzyme", "enzymes", "enzymatic"],
	("biology*", "molecular"): ["molecular", "biomolecular"],
	("biology*", "microbiology"): ["microbiology", "microbiologist"],
	("biology*", "stem cells"): ["stem cells", "stem cell"],
	("biology*", "organ systems"): ["organ systems", "organ system"],
	("biology*", "neuropeptides"): ["neuropeptides", "neuropeptide"],
	("biology*", "proteins"): ["protein", "proteins", "protein-protein", "proteomics"],
	("biology*", "RNA"): ["RNA", "ribosomal", "ribosome", "RNA-protein"],
	("biology*", "DNA"): ["DNA"],
	("biology*", "CRISPR"): ["CRISPR-Cas9", "CRISPR", "CRISPR-Cas"],
	("biology*", "cell death"): ["cell death", "ferroptosis"],
	("biology*", "cell"): ["cell", "cells", "cellular", "cell-based", "single-cell"],
	("biology*", "eukaryotes"): ["eukaryotes", "eukaryotic"],
	("biology*", "pathobiology"): ["pathobiology"],
	("biology*", "microbiome"): ["microbiome", "microbe", "microbes", "microbial", "host-microbe", "microbiota"],
	("biology*", "gut"): ["gut", "gastric", "gastroenterology", "gastrointestinal", "GI"],
	("biology*", "pathophysiology"): ["pathophysiology"],
	("biology*", "paleobiology"): ["paleobiology"],

	# Ecology
	"ecology*": ["ecology", "ecological"],
	("ecology*", "microecology"): ["microecology"],
	("ecology*", "macroecology"): ["macroecology"],

	# Genetics
	("genes*", "genetics"): ["genetics", "genetic"],
	("genes*", "DNA"): ["DNA"],  
	("genes*", "epigenetics"): ["epigenetics", "epigenetic", "epigenome"], 
	("genes*", "genomics"): ["genomics", "genomic", "genome", "immunogenomics", "pharmacogenomics"], 
	("genes*", "gene delivery"): ["gene delivery"], 
	("genes*", "gene expression"): ["gene expression"],
	("genes*", "proteomics")

	# Pharmacology
	"pharmacology*": ["pharmaceutical", "pharmacogenetics", "pharmacological", "pharmacology"], 
	("pharmacology*", "neuropharmacology"): ["neuropharmacology"],
	("pharmacology*", "antidepressants"): ["antidepressants","antidepressant"],
	("pharmacology*", "drugs"): ["drugs", "drug"],
	("pharmacology*", "PKPD"): ["PKPD", "Pharmacokinetic-Pharmacodynamic"],
	("pharmacology*", "QSP"): ["QSP", "Quantitative Systems Pharmacology"],
	("pharmacology*", "PBPK"): ["PBPK", "Physiologically based pharmacokinetic"],
	("pharmacology*", "pharmacogenomics"): ["pharmacogenomics"],
	("pharmacology*", "pharmacogenetics"): ["pharmacogenetics"],
	("pharmacology*", "drug discovery"): ["drug discovery"],
	("pharmacology*", "drug development"): ["drug development"],
	
	# Cancer
	("cancer*", "cancer"): ["cancer", "cancerous", "carcinogenesis"], 
	("cancer*", "oncology"): ["oncology", "oncologist"],
	("cancer*", "radiotherapy"): ["radiotherapy", "radiation therapy"],
	("cancer*", "tumor"): ["melanoma", "melanomas", "tumor", "tumors", "pro-tumorigenic", "tumorigenic"],

	# More Health
	"cardiology*": ["cardiology", "cardiologist", "cardiac surgery", "cardiovascular", "heart"],
	"sexual health*": ["sexually", "STI", "STIs", "sexual"], 
	"medical images*": ["medical images", "medical imaging"],
	("medical images*", "MRI"): ["MRI"],
	("medical images*", "fMRI"): ["fMRI"],

	# Marine
	"marine*": ["marine", "sea", "ocean", "oceans"],
	("marine*", "sea turtles"): ["sea turtle", "sea turtles"],
	("marine*", "deep-sea"): ["deep-sea", "deep sea"],
	("marine*", "oceanography"): ["oceanography", "oceanographic", "oceanographer"],

	# Fish
	"fish*": ["fish", "fishes"],
	("fish*", "fisheries"): ["fisheries", "fishery"],

	# Animals
	"animals*": ["animals", "animal", "wildlife"],
	("animals*", "veterinary"): ["vet", "veterinarian", "veternary"],
	("animals*", "zoology"): ["zoology", "zoological"], 
	("animals*", "primates"): ["primates", "primate", "chimpanzee", "chimpanzees", "bonobo", "bonobos"], 
	
	# Chemistry
	"chemistry*": ["chemistry", "chemist"],
	("chemistry*", "protein chemistry"): ["protein chemistry"],
	("chemistry*", "physical chemistry"): ["physical chemistry"],    
	("chemistry*", "catalysts"): ["catalysts", "catalyst", "catalysis"],
	("chemistry*", "toxicology"): ["toxicology", "toxicological"],
	("chemistry*", "radiochemistry"): ["radiochemistry", "radiochemical"],
	("chemistry*", "chemometrics"): ["chemometrics"],
	("chemistry*", "mechanochemistry"): ["mechanochemistry"],  
	("chemistry*", "electrochemistry"): ["electrochemistry", "electrochemical"], 
	("chemistry*", "nanochemistry"): ["nanochemistry", 'nanochemicals'],  
	
	# Geoscience
	"geoscience*": ["geology", "geoscience", "geosciences", "geoscientist"],
	("geoscience*", "geochemistry"): ["geochemistry", "geochemist"],
	("geoscience*", "biogeochemistry"): ["biogeochemistry", "biogeochemical"],
	("geoscience*", "geophysics"): ["geophysics", "geophysicist", "geophysical"],  
	("geoscience*", "mineralogy"): ["mineralogy", "minerologist"],
	("geoscience*", "geobiology"): ["geobiology", "geobiologist"], 

	# Eco-friendly
	("eco-friendly", "sustainability"): ["sustainability", "sustainable"], 
	("eco-friendly", "conservation"): ["conservation", "conservational", "conservationist"],
	("eco-friendly", "preservation"): ["presservation", "preservational"],
	("eco-friendly", "renewable energy"): ["renewable energy", "renewable energies"],

	# Mathematics
	"mathemtics*": ["mathematics", "mathematical", "math"],
	("mathemtics*", "complex differential"): ["complex differential"],
	("mathemtics*", "graph theory"): ["graph theory"], 
	("mathemtics*", "number theory"): ["number theory"], 
	("mathemtics*", "algebra*"): ["algebra", "algabreic"],
	("mathemtics*", "geometry*"): ["gemeomtry", "gemoetric"],
	("mathemtics*", "combinatorics"): ["combinatorics"],
	("mathemtics*", "math modeling"): ["math modeling", "mathematical modeling"],

	# Computer Science
	"computer science*": ["computer science", "CS", "computer scientist", "computing"],
	("computer science*", "supercomputers"): ["supercomputers"], 
	("computer science*", "complexity theory"): ["complexity theory", "computational complexity"],
	("computer science*", "compilers"): ["compilers", "compiler"],
	("computer science*", "high-performance computing"): ["high-performance computing"], 
	("computer science*", "affective computing"): ["affective computing", "emotion recognition"],
	("computer science*", "cryptography"): ["cryptography"],
	("computer science*", "social computing"): ["social computing"],
	("computer science*", "scientific computing"): ["scientific computing"],
	("computer science*", "coding"): ["coding", "programming"],
	("computer science*", "computer vision"): ["computer vision"],
	("computer science*", "computer architecture"): ["computer architecture"],
	("computer science*", "computer networks"): ["computer networks"],

	# Data Science
	"data science*": ["data science", "data scientist"],
	("data science*", "data mining") : ["data mining", "text mining"],
	("data science*", "data analysis") : ["data analysis"],
	("data science*", "statistics"): ["statistics", "stats"],
	("data science*", "databases"): ["databases", "databases"],
	("data science*", "informatics"): ["informatics"],
	("data science*", "information theory"): ["information theory"],


	# Signal Processing
	"signal processing*": ["signal processing"],
	("signal processing*", "PCA"): ["principal component analysis", "PCA"],
	("signal processing*", "feature extraction"): ["feature extraction"],
	("signal processing*", "image processing"): ["image processing"],

	# Physics
	"physics*": ["physics", "physicist"],
	("physics*", "optoelectronics"): ["optoelectronics", "optoelectronic"], 

	# Quantum
	"quantum*": ["quantum"],
	("quantum*", "quantum gravity"): ["quantum gravity"], 
	("quantum*", "quantum field"): ["quantum field"],
	("quantum*", "quantum mechanics"): ["quantum mechanics", "quantum mechanical"],  
	("quantum*", "quantum-chemical"): ["quantum-chemical"], 
	("quantum*", "post-quantum cryptography"): ["post-quantum cryptography"], 

	# Space (Physics)
	"space*": ["outer space", "space", "astro"],
	("space*", "planets"): ["planets", "planet", "planetary"],
	("space*", "exoplanet"): ["exoplanet", "exoplanets"],  
	("space*", "astronomy"): ["astronomy"],  
  	("space*", "astrobiology"): ["astrobiology", "exobiology"], 
  	("space*", "telescope"): ["telescope", "telescopes"],  
  	("space*", "stars"): ["star", "stars"],
  	("space*", "galaxy"): ["galaxy", "galaxies"],
  	("space*", "spaceflight"): ["spaceflight"],
  	("space*", "dark matter"): ["dark matter"],
  	("space*", "neutron stars"): ["neutron stars", "neutron star"],
  	("space*", "X-ray binary stars"): ["X-ray binary stars", "X-ray binary star"],
  	("space*", "astrophysics"): ["astrophysics", "astrophysicist", "astrophysical"],
  	("space*", "cosmology"): ["cosmology", "cosmological"],

  	# Aerospace
  	"aero*": ["aerospace", "aeronautical"],
  	("aero*", "aerodynamics"): ["aerodynamics", "aerodynamic"],

	# Archeology
	"archeology*": ["archaeological", "archeology"],
	("archeology*", "bioarchaeology"):["bioarchaeology"],
	("archeology*", "archaeobotanical"):["archaeobotanical"],

	# Anthropology
	"anthropology*": ["anthropology", "anthropologist"], 
	("anthropology*", "biological anthropology"): ["biological anthropology"],
	("anthropology*", "forensic anthropology"): ["forensic anthropology"],
	("anthropology*", "cultural anthropolgy"): ["cultural anthropology"],

	# AI / ML
	("AI/ML*", "AI"): ["AI", "artificial intelligence"],
	("AI/ML*", "AGI"): ["AGI", "artificial general intelligence"],
	("AI/ML*", "ML"): ["ML", "machine learning"],
	("AI/ML*", "deep learning"): ["deep learning"],
	("AI/ML*", "neural networks"): ["neural networks", "neural nets"],

	# AR/VR
	("AR/VR*", "AR"): ["augmented reality", "AR"],
	("AR/VR*", "VR"): ["virutal reality", "VR"],
	("AR/VR*", "XR"): ["extended reality", "XR"],

	# Nanotechnology 
	("nano*", "nanotechnology"): ["nanotechnology", "nanotech"],
	("nano*", "nanoparticles"): ["nanoparticles", "nano-particles"],
	("nano*", "nanoscience"): ["nanoscience"],
	("nano*", "nanophotonics"): ["nanophotonics", "nanophotonic"],
	("nano*", "nanomaterials"): ["nanomaterials", "nanomaterial"],
	("nano*", "matierals science"): ["matierals science"], 

	# Plurals and Other Variations
	"technology": ["technology", "tech", "technologies"],
	"biochemistry": ["biochemistry", "biochemical", "biochem", "chemical biology"], 
	"psychology": ["psychology", "psychological", "biopsychology"],
	"anesthesiology": ["anesthesiology", "anesthesiologist"],
	"immunity*": ["immunology", "immunogenomics", "immune", "immune-mediated", "autoimmune", "autoimmunity",
	"neuroimmunology", "immunotherapy", "immune-mediated", "immunity"],
	"infection": ["infection", "infections", "infectuous"],
	"robotics": ["robot", "robots", "robotics", "robotic"],
	"diagnosis": ["diagnosis", "diagnoses", "diagnose", "diagnostic"],
	"medical device": ["medical devices", "medical device"],
	"bioinformatics": ["bioinformatics", "biomedical informatics", "biostatistics"],
	"bacteria*": ["bacteria", "bacterium", "bacterial"], 
	"vaccines*": ["vaccines", "vaccine", "vaccinations", "vaccination", "gene-vaccine"], 
	"infectious*": ["infectious", "infection", "infections"], 
	"intervention": ["intervention", "interventional"],
	"antimicrobial": ["antimicrobials", "antimicrobial"],
	"antibiotics": ["antibiotic", "antibiotics", "antibiotic-resistant"],
	"sequencing": ["sequence", "sequencing"],
	"HCI*": ["human-machine interaction", "human-computer interaction"], 
	"BMI*": ["brain-machine interface", "brain-machine interfaces", "brain-computer interface", "brain-computer interfaces", 
			"BMI", "BMIs", "brain-computer interfaces", "brain-computer interface", "neural interfacing", "neural interface"], 
	"NLP": ["NLP", "natural language processing"],
	"data engineering": ["data engineering", "data engineer"], 
	"statistics": ["statistics", "stats", "statistical"],
	"medicine*": ["medicine", "medicinal", "medical", "health", "healthcare", "medical sciences", "pre-med"],  
	"plants*": ["plants", "plant", "botony", "botanical"],
	"exercise*": ["fitness", "exercise science"], 
	"prosthetics*": ["prosthetics", "prostheses", "bionic"], 
	"environment": ["environmental", "environment"],
	"biomarkers": ["biomarkera", "biomarkers", "biomarker", "bio markers"],
	"parasites": ["parasite", "parasites", "parasitic"], 
	"nutrition": ["nutritional science", "nutritional", "nutrition"], 
	"rehabilitation": ["rehabilitation", "rehabilitative", "rehab"],
	"security*": ["security", "cybersecurity", "data privacy", "privacy-preserving", "privacy"],
	"biofilm": ["biofilm", "biofilms"],
	"atoms": ["atoms", "atom", "atomic"],
	"developmental": ["developmental", "development"],
	"treatment*": ["treat", "treatment", "treatments", "therapy", "therapies", "therapeutic", "therapeutics"],
	"ultrasounds":["ultrasound", "ultrasounds"],
	"algorithm": ["algorithm", "algorithms"],
	"computer": ["computer", "computers"],
	"evolution": ["evolution", "evolutionary"], 
	"neutrons": ["neutrons", "neutron"],
	"electrons": ["electrons", "electron"],
	"protons": ["protons", "proton"],
	"isotopes": ["isotopes", "isotope"], 
	"polymers": ["polymers", "polymer"],
	"carbon": ["carbon", "carbon-carbon"],
	"OBGYN": ["obsectrics", "obstetrician", "gynecology", "gynecologist"],
	"strokes": ["stroke", "strokes"],
	"amino acids": ["amino acids", "amino acid"],
	"image recognition": ["image recognition", "facial recognition"],
	"classification": ["classification", "classifier", "classify", "classifying"],
	"analytics": ["analytics", "analytic", "analytical"],
	"metabolism": ["metabolism", "metabolic"],
	"hurricanes": ["hurricanes", "hurricane"],
	"environmental science": ["environmental science", "environmental sciences"], 
	"system": ["system", "systems"],
	"collision": ["collisions", "collision"],
	"retina": ["retina", "retinal"],
	"tenure": ["tenure", "tenured"],
	"biotechnology": ["biotech","biotechnology"],
	"minerals": ["minerals", "mineral"],
	"neurons": ["neurons", "neuron", "neuronal", "glia", "glial"],
	"orthopedics": ["orthopedics", "orthopedic", "orthopaedics", "orthopaedic"],
	"serotonin": ["serotonin", "serotonergic"],
	"hydrodynamics": ["hydrodynamics", "hydrodynamic"],
	"technology": ["technology", "tech"],
	"photovoltaics": ["photovoltaics", "photovoltaic", "PV"],
	"exosuits": ["exosuits", "exosuit"],
	"drones": ["drones", "drone"],
	"low-dimensional": ["low-dimensional", "low dimensional"],
	"nursing": ["nursing", "nurse", "nurses"],
	"pediatric": ["pediatrics", "pediatric"],
	"phosphatase": ["phosphatase", "phosphatases"],
	"forensics": ["forensics", "forensic"],
	"alloys": ["alloys", "alloy"],
	"high-temperature": ["high-temperature", "high temperature"],
	"organisms": ["organisms", "organism", "organismal", "microorganism", "microorganisms"],
	"wildfires": ["wildfires", "wildfire"],
	"FPOP": ["FPOP", "Fast photochemical oxidation of proteins"],
	"macrophages": ["macrophages", "macrophage"],
	"glaciers": ["glacier", "glaciers"]
}

hobbies = ["sports", "board games", "piano", "guitar", "violin", "cello", "flute", "saxaphone",
	"viola", "trumpet", "clarinet", "trumbone", "drums", "ukulele", "cooking",
	"baking", "outreach", "snowboarding", "surfing", "mountain biking", "science fiction", "climbing", 
	"squash", "basketball", "Netflix", "baseball", "softball", "hockey", "field hockey", "ice hockey", 
	"figure skating", "documentaries", "camping", "ballet", "jazz", "singing", "salsa", "bachata",
	"foreign languages", "yoga", "traveling", "pilates", "running", "gymnastics","backpacking", 
	"rock climbing", "bouldering", "acting", "theatre", "concerts", "female empowerment", "skiing", "drawing", 
	"painting", "kayaking", "canoing", "knitting", "skateboarding", "ice skating", "snowboarding",
	"audiobooks", "podcasts", "tv", "movies", "swimming", "music", "rafting", "spanish", "german", "french", 
	"romanian", "arabic", "photography", "gardening", "soccer", "coding", "racquet", "food", "science fiction", 
	"chess", "fishing", "tennis", "football", "golf", "cycling", "rugby", "netball", "pottery", "sculpting", 
	"badminton", "anime", "public service", "ethics", "economics", "policy", "diplomacy", "arms control"] 

hobbies_dict = {"reading": ["reading", "read"], "writing": ["writing", "write"], "dogs": ["dogs", "dog"], 
	"cats": ["cats", "cat"], "video games": ["gaming", "video games"], "volunteering": ["volunteering", "volunteer"],
	"dancing": ["dancing", "dance"], "equestrian": ["horseback riding", "equestrian"], "crew": ["rowing", "crew"],
	"hip hop": ["hip-hop", "hip hop"], "hiking": ["hiking", "hike"], "art": ["artwork", "artworks", "art", "arts"]
}

fellowships = ["Marschall", "Churchill", "Rhodes", "Fulbright", "Soros", "OxCam", "Truman", "Boren", "Udall"]

fellowships_dict = {
	"fellowship": ["fellowship", "fellowships", "fellow", "scholarship"], 
	"NSF GRFP": ["NSF", "GRFP", "National Science Foundation Graduate"],
	"NDSEG": ["NDSEG", "National Defense Science Engineering Graduate"], 
	"Gates": ["Gates-Cambridge", "Gates", "Gates Cambridge"], 
	"NIH F31": ["NIH F31", "National Research Service Award"]}

career = ["company", "CDC", "NIH", "NASA", "national", "NOAA", "USDA", "academia", "professor", "MD/PhD", "research fellow", "industry"]

career_dict = {
	"lab": ["lab", "laboratory"], 
	"government": ["government", "governmental"], 
	"post-doc": ["postdoc", "post-doc","post-doc", "postdoctoral"],
	"startups": ["startups", "startup"],
	"consulting": ["consulting", "consultant", "consultancy"],
	"entrepreneur": ["entrepreneurial", "entrepreneur"],
	"visiting assistant": ["visiting assistant", "VAP"],
	"management": ["manager", "management"],
	"masters": ["M.S.", "Masters"]
	"MPhil": ["MPhil", "Master of Philosophy"]
	"PhD": ["PhD", "Ph.D."]
}

primary_fields = ["Chemistry", "Computer & Information Sciences & Engineering", "Engineering", "Geosciences", 
	"Life Sciences", "Materials Research", "Mathematical Sciences", "Medicine", "Physics & Astronomy", "Psychology"]

sub_fields = [

	# Biology
	"Chemical Catalysis", "Macromolecular, Supramolecular, and Nanochemistry", "Chemical Measurement and Imaging", 
	"Chemical Structure, Dynamics, and Mechanism", "Chemical Theory, Models, and Computational Methods", 
	"Environmental Chemical Systems", "Sustainable Chemistry", "Chemistry of Life Processes", "Chemical Synthesis", 

	# Computer & Information Sciences & Engineering
	"Bioinformatics and Other Informatics", "Data Mining and Information Retrieval Databases", "Graphics and Visualization", 
	"Human-Computer Interaction", "Machine Learning", "Natural Language Processing", "Robotics and Computer Vision", 
	"Algorithms and Theoretical Foundations", "Communication and Information Theory", "Computational Science and Engineering", 
	"Computer Architecture", "Computer Networks", "Computer Security and Privacy", "Computer Systems and Embedded Systems", 
	"Formal Methods, Verification, and Programming Languages", "Software Engineering", 

	# Engineering 
	"Aeronautical and Aerospace Engineering", "Energy Engineering", "Nuclear Engineering", "Optical Engineering", "Systems Engineering", 
	"Bioengineering", "Biomedical Engineering", "Chemical Engineering", "Polymer Engineering", "Civil Engineering", "Environmental Engineering", 
	"Ocean Engineering", "Computer Engineering", "Electrical and Electronic Engineering", "Industrial Engineering and Operations Research", 
	"Materials Engineering", "Mechanical Engineering",

	# Geosciences
	"Aeronomy", "Atmospheric Chemistry", "Climate and Large-Scale Atmospheric Dynamics", "Magnetospheric Physics", "Paleoclimate", 
	"Physical and Dynamic Meteorology", "Solar Physics", "Geobiology", "Geochemistry", "Geodynamics", "Geomorphology", "Geophysics", 
	"Glaciology", "Hydrology", "Paleontology and Paleobiology", "Petrology", "Sedimentary Geology", "Tectonics", "Biogeochemistry", 
	"Biological Oceanography", "Chemical Oceanography", "Marine Biology", "Marine Geology and Geophysics", "Physical Oceanography", 

	# Life Sciences
	"Biochemistry", "Biophysics", "Structural Biology", "Cell Biology", "Ecology", "Environmental Biology", "Biodiversity", 
	"Evolutionary Biology", "Systematics", "Bioinformatics and Computational Biology", "Genetics", "Genomics", "Proteomics", 
	"Microbial Biology", "Molecular Biology", "Systems Biology", "Neurosciences", "Developmental Biology", "Organismal Biology", 
	"Physiology", 

	# Materials Research
	"Biomaterials", "Ceramics", "Chemistry of Materials", "Electronic Materials", "Materials Theory", "Metallic Materials", 
	"Photonic Materials", "Physics of Materials", "Polymers", 

	# Mathematical Sciences
	"Algebra, Number Theory, and Combinatorics", "Analysis", "Geometric Analysis", "Logic or Foundations of Mathematics", 
	"Probability", "Statistics", "Topology", "Applied Mathematics", "Biostatistics", "Computational and Data-enabled Science",
	"Computational Mathematics", "Computational Statistics", "Mathematical Biology", 

	# Physics & Astronomy
	"Astronomy and Astrophysics", "Atomic, Molecular, and Optical Physics", "Nuclear", "Plasma", "Condensed Matter Physics", 
	"Particle Physics", "Physics of Living Systems", "Solid State", "Theoretical Physics", 

	# Psychology
	"Cognitive Psychology", "Cognitive Neuroscience", "Computational Psychology", "Psycholinguistics", "Developmental", 
	"Experimental or Comparative", "Neuropsychology", "Perception and Psychophysics", "Physiological","Quantitative"]

# Keyword Processor for Mentor and Mentee Profile Data
keyword_processor = KeywordProcessor()

# Adds keywords from dictionaries
keyword_processor.add_keywords_from_dict(science_dict)
keyword_processor.add_keywords_from_dict(hobbies_dict)
keyword_processor.add_keywords_from_dict(fellowships_dict)
keyword_processor.add_keywords_from_dict(career_dict)

# Adds keywords from lists
keyword_processor.add_keywords_from_list(science)
keyword_processor.add_keywords_from_list(hobbies)
keyword_processor.add_keywords_from_list(fellowships)
keyword_processor.add_keywords_from_list(career)
keyword_processor.add_keywords_from_list(primary_fields)
keyword_processor.add_keywords_from_list(sub_fields)


# Experimenting

# print(keyword_processor_mentee.get_all_keywords())

mentee_keys = keyword_processor.extract_keywords(" physical chemistry Hi all, I'm neural networks Quinn and I am a current junior at Bowling Green State University in Ohio. I study Biochemistry but my research is in bioinformatics where I work on RNA-protein sequences and interaction prediction from three-dimensional structures! I plan on pursuing a Ph.D. in Computational Biology to continue this type of sequence work, although I am certainly open to new topics of research. If you are interested in bioinformatics I would love to get to know you! There is not a bioinformatics program at my school so it is rare that I get to talk to other bioinformaticians in the making ")

mentor1 = keyword_processor.extract_keywords("Hi, everyone! computational oncology quantum mechanical nanoparticles nanotechnology I’m Cody, a third year biochemistry and genetics double major with double minors in bioinformatics and statistics (always a mouthful) at Texas A&M University. I’m part of the Beckman Scholars research program, and for my research, I study mechanisms of bacteriophage lysis, which is how bacterial viruses escape from the host cell after infection. Nice to meet all of you!!")

mentor2 = keyword_processor.extract_keywords("Hi y'all! neuropharmacology I'm Noah, and I'm a sophomore studying physics and math at North Carolina State University. I plan to pursue a PhD in astrophysics, with a focus on theoretical and/or computational astrophysics-- right now, I'm having a great time blowing up stars on computers to push our understanding of fundamental physics. Super excited to meet everyone!!")

mentor3 = keyword_processor.extract_keywords("Hi all! : cultural anthropology: I’m Vaishnavi, a sophomore at MIT studying Computer Science & Molecular Biology. My research is in basic cancer cell biology; I work on a recently discovered type of cell death called ferroptosis that can potentially be used to eliminate aggressive cancer cells. I’m hoping to pursue an MD/PhD in cancer and/or cell biology - excited to get to know you all!")

# print("Mentor 1 ", mentor1)
# print("Mentor 2 ", mentor2)
# print("Mentor 3 ", mentor3)

# Transforms list of mentee keys so that tuple keys (general_key, specific_key) are split into two separate strings
# Returns final list without duplicates 
def general_and_specific(keys):
	mentee = []
	for key in keys:
		if type(key) == tuple:
			key1, key2 = key
			mentee.append(key1)
			mentee.append(key2)
		else:
			mentee.append(key)
	return list(set(mentee))

# Returns number of common keys between two lists of keys 
def common_key(a, b):
	num_common = 0
	a_set = set(a)
	b_set = set(b)
	if (a_set & b_set):
		num_common = len(a_set & b_set)
	return num_common

mentee_keys = general_and_specific(mentee_keys)
print("\n", mentee_keys)
print("\n", general_and_specific(mentor1))
print("\n", general_and_specific(mentor2))
print("\n", general_and_specific(mentor3))

print(common_key(mentee_keys, general_and_specific(mentor1)))
print(common_key(mentee_keys, general_and_specific(mentor2)))
print(common_key(mentee_keys, general_and_specific(mentor3)))


# How It Will Work:
# have huge keyword list 
# for each mentee
	# get list of keywords
	# make sure these two have not been paired before
		# if mentor not in mentee.previous_mentors 
	# compare to all potential mentors that have the specific advice category that the mentee selected
		# if advice in mentor.advise_areas
			# use the mentee keyword list to search for matching keywords in mentors
# take out mentors from list when they've found a match so they aren't used in future searches (for as many matches as they're willing)

# Things to Consider:
# 1. Age/career level group that a mentee wants to match with
	# junior, senior, grad, post-doc, career/professional, etc
# 2. Advise areas that a mentor is willing to work in
# 3. Advise area that the mentee wants
# 4. Number of mentees a mentor is willing to take on (provide them with a list of that many top matches)
# 5. Previous mentor matches

# Profile Instructions
# Allow NO slashes (/) s in writing! 
# Write out full instutution names i.e. "Massachus"

# Field of study, research area, career goals/grad school goals, instituions, 




