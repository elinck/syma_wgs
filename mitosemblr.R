#### mitosemblR: an R wrapper for automated mitochondrial genome assembly from off-target reads ###

blastR <- function(read_dir, reference){
  R1 <- list.files(read_dir, full.names=T) %>% grep("R1_matches",.,value=T)
  R2 <- list.files(read_dir, full.names=T) %>% grep("R2_matches",.,value=T)
  ST <- list.files(read_dir, full.names=T) %>% grep("singleton",.,value=T)
  commands <- c()
  for(i in 1:20){
    system('echo "$@"_starting blatq search to pull seeds for PRICE')
    commands[i] <- paste0("blatq -t=dna", reference, R1[i], 
                          R1[i] %>% basename() %>% strsplit("_L006_R._001") %>% unlist() %>% .[1],"R1_matches.m8;", 
                          "blatq -t=dna", reference, R2[i], 
                          R2[i] %>% basename() %>% strsplit("_L006_R._001") %>% unlist() %>% .[1],"R2_matches.m8;", 
                          "blatq -t=dna", reference, ST[i],
                          ST[i] %>% basename() %>% strsplit("_L006_R._001") %>% unlist() %>% .[1],"singleton_matches.m8")
    system('echo "$@" Blatq done')
  }
  for(i in commands){
    system(i)
  }
}

excerptR <- function(read_dir, m8_dir, excerptByIDPath){
  R1 <- list.files(read_dir, full.names=T) %>% grep("R1_001",.,value=T)
  R2 <- list.files(read_dir, full.names=T) %>% grep("R1_002",.,value=T)
  ST <- list.files(read_dir, full.names=T) %>% grep("singleton",.,value=T)
  M1 <- list.files(m8_dir, full.names=T) %>% grep("R1_matches",.,value=T)
  M2 <- list.files(m8_dir, full.names=T) %>% grep("R2_matches",.,value=T)
  MS <- list.files(m8_dir, full.names=T) %>% grep("singleton_matches",.,value=T)
  commands <- c()
  for(i in 1:20){
    system('echo "$@"Now starting to excerpt reads with blatq hits by ID')
    commands[i] <- paste0(excerptByIDPath, R1[i], M1[i], ">", 
                          R1[i] %>% basename() %>% strsplit("R1_matches.m8") %>% unlist() %>% .[1],"R1_matches.fq;",
                          excerptByIDPath, R2[i], M2[i], ">", 
                          R2[i] %>% basename() %>% strsplit("R2_matches.m8") %>% unlist() %>% .[1],"R2_matches.fq;",
                          excerptByIDPath, ST[i], MS[i], ">", 
                          ST[i] %>% basename() %>% strsplit("singleton_matches.m8") %>% unlist() %>% .[1],"singleton_matches.fq;")
  }
    system('"$@"done Excerpting reads...')
  for(i in commands){
    system(i)
  }
}

seedR <- function(read_dir){
  R1 <- list.files(read_dir, full.names=T) %>% grep("R1_matches.fq",.,value=T)
  R2 <- list.files(read_dir, full.names=T) %>% grep("R2_matches.fq",.,value=T)
  ST <- list.files(read_dir, full.names=T) %>% grep("singleton_matches.fq",.,value=T)
  for(i in 1:20){
    commands[i] <- paste0('cat', R1[i], R2[i], ST[i], '>', 
                          R1[i] %>% basename() %>% strsplit("R1_matches.m8") %>% unlist() %>% .[1],'SEEDS.fastq')
  }
  for(i in commands){
    system(i)
  }
}

SPAdeR <- function(read_dir, spades_path){
  seeds <- list.files(read_dir, full.names=T) %>% grep("SEEDS.fastq",.,value=T)
  system('echo "$@"Now starting SPAdes assembly...')
  for(i in 1:20){
    commands[i] <- paste0('python', spades_path, '--careful -s', seeds[i], 
                          seeds[i] %>% basename() %>% strsplit("SEEDS.m8") %>% unlist() %>% .[1],'SEEDS_assembled')
  }
}


priceR <- function(read_dir, price_path, insertSize, identityPercentage, cycles, mol, mpi, MPI, logpath){
  R1 <- list.files(read_dir, full.names=T) %>% grep("R1_matches.fq",.,value=T)
  R2 <- list.files(read_dir, full.names=T) %>% grep("R2_matches.fq",.,value=T)
  ST <- list.files(read_dir, full.names=T) %>% grep("singleton_matches.fq",.,value=T)
  seeds <- list.files(read_dir, full.names=T) %>% grep("SEEDS_assembled",.,value=T)
  system('echo "$@"Now starting PRICE assembly using reads')
  for(i in 1:20){
    commands[i] <- paste(price_path, '-fpp', R1[i], R2[i], insertSize, identityPercentage,
                         '-icf', seeds, '1 1 5', '-nc', cycles, '-mol', mol, '-mpi', mpi,
                         '-MPI', MPI, '-a', '-a 8 -target 85 1 1 1 -o', 
                         seeds[i] %>% basename() %>% strsplit("SEEDS_assembled") %>% unlist() %>% .[1],'mtDNA.fasta', '-o',
                         seeds[i] %>% basename() %>% strsplit("SEEDS_assembled") %>% unlist() %>% .[1],'mtDNA.fasta',
                         '-maxHp 25 -logf', logpath)
  }
}
