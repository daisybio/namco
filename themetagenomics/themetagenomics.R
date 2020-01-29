#credit to themetagenomics vis_taxa.R
prepare_vis <- function(object,taxa_bar_n=30,top_n=7,method='huge',corr_thresh=.01,lambda_step=.01,...){
  
  topics <- object$topics
  topic_effects <- object$topic_effects
  tax_table <- object$topics$tax_table
  
  method <- match.arg(method)
  
  fit <- topics$fit
  K <- fit$settings$dim$K
  vocab <- fit$vocab
  taxa <- tax_table[vocab,]
  taxon <- paste0(pretty_taxa_names(taxa),' (',vocab,')')
  taxa <- rename_taxa_to_other(topics$docs,taxa,top_n=top_n,type='docs',as_factor=TRUE)
  rownames(taxa) <- taxon
  
  corr <- stm::topicCorr(fit,method=method,cutoff=corr_thresh)
  
  colors_order <- sort(as.character(c(-1,0,1)))
  colors <- c('gray','indianred3','dodgerblue3','indianred4','dodgerblue4','gray15')
  names(colors) <- c('0','1','-1','2','-2','00')
  
  covariates <- lapply(names(topic_effects),identity)
  names(covariates) <- tolower(names(topic_effects))
  est_range <- range(c(sapply(topic_effects, function(x) unlist(x$est))))
  
  doc_lengths <- sapply(topics$docs,function(x) sum(x[2,]))
  
  theta <- fit$theta
  rownames(theta) <- names(topics$docs)
  colnames(theta) <- paste0('T',1:K)
  if (min(theta) == 0){
    min_theta <- min(theta[theta > 0])
    theta <- theta + min_theta
    theta <- theta/rowSums(theta)
  }
  
  beta <- exp(fit$beta$logbeta[[1]])
  colnames(beta) <- taxon
  rownames(beta) <- paste0('T',1:K)
  if (min(beta) == 0){
    min_beta <- min(beta[beta > 0])
    beta <- beta + min_beta
    beta <- beta/rowSums(beta)
  }
  
  topic_freq <- colSums(theta*doc_lengths)
  topic_marg <- topic_freq/sum(topic_freq)
  term_topic_freq <- beta*topic_freq
  
  # compute term frequencies as column sums of term.topic.frequency
  # we actually won't use the user-supplied term.frequency vector.
  # the term frequencies won't match the user-supplied frequencies exactly
  # this is a work-around to solve the bug described in Issue #32 on github:
  # https://github.com/cpsievert/LDAvis/issues/32
  term_freq <- colSums(term_topic_freq) # topic_freq <- colSums(theta*doc_lengths)
  
  # marginal distribution over terms (width of blue bars)
  term_marg <- term_freq/sum(term_freq)
  
  beta <- t(beta)
  
  # compute the distinctiveness and saliency of the terms:
  # this determines the R terms that are displayed when no topic is selected
  topic_g_term <- beta/rowSums(beta)  # (W x K)
  kernel <- topic_g_term*log(t(t(topic_g_term)/topic_marg))
  distinctiveness <- rowSums(kernel)
  saliency <- term_marg*distinctiveness
  
  # Order the terms for the "default" view by decreasing saliency:
  default_terms <- taxon[order(saliency,decreasing=TRUE)][1:taxa_bar_n]
  counts <- as.integer(term_freq[match(default_terms,taxon)])
  taxa_n_rev <- rev(seq_len(taxa_bar_n))
  
  default <- data.frame(Term=default_terms,
                        logprob=taxa_n_rev,
                        loglift=taxa_n_rev,
                        Freq=counts,
                        Total=counts,
                        Category='Default',
                        stringsAsFactors=FALSE)
  default$Taxon <- taxa[default$Term,'Phylum']
  default$Term <- factor(default$Term,levels=rev(default$Term),ordered=TRUE)
  
  topic_seq <- rep(seq_len(K),each=taxa_bar_n)
  category <- paste0('T',topic_seq)
  lift <- beta/term_marg
  
  # Collect taxa_bar_n most relevant terms for each topic/lambda combination
  # Note that relevance is re-computed in the browser, so we only need
  # to send each possible term/topic combination to the browser
  find_relevance <- function(i){
    relevance <- i*log(beta) + (1 - i)*log(lift)
    idx <- apply(relevance,2,function(x) order(x,decreasing=TRUE)[seq_len(taxa_bar_n)])
    indices <- cbind(c(idx),topic_seq)
    data.frame(Term=taxon[idx],
               Category=category,
               logprob=round(log(beta[indices]),4),
               loglift=round(log(lift[indices]),4),
               stringsAsFactors=FALSE)
  }
  
  lambda_seq <- seq(0,1,by=lambda_step) # 0.01
  tinfo <- lapply(as.list(lambda_seq),find_relevance)
  tinfo <- unique(do.call('rbind', tinfo))
  tinfo$Total <- term_freq[match(tinfo$Term,taxon)]
  rownames(term_topic_freq) <- paste0('T',seq_len(K))
  colnames(term_topic_freq) <- taxon
  tinfo$Freq <- term_topic_freq[as.matrix(tinfo[c('Category','Term')])]
  tinfo <- rbind(default[,-7],tinfo)
  
  new_order <- default
  
  return (out = list(beta = beta, 
                     tinfo = tinfo,
                     taxa_bar_n = taxa_bar_n,
                     taxa = taxa,
                     default = default,
                     topic_marg = topic_marg,
                     corr = corr,
                     lambda_step = lambda_step,
                     covariates = covariates,
                     colors=colors))
}