#' Estimate ICCs, correlations, and descriptive Statistics
#'
#' This function creates the ICC matrix
#' @param x Data object.
#' @param group Nesting variable.
#' @param title Table caption.
#' @param gmc Provide group-mean center correlations.  Default is FALSE.
#' @param alpha.order Sort variables alphabetically.  Default is False
#' @param result Output options.  Default is html table to viewer.  Option "text" returns output to console only.
#' @return A correlation table with sample statistics and ICC estimates
#' @export

icc.corrs <- function(x, group, title = "Descriptive Stats", gmc = FALSE, alpha.order = FALSE, result = "html") {

  # ungroup if original dataset is grouped
  if (dplyr::is.grouped_df(x)) {

    x <- x %>% dplyr::ungroup()
  }

  # rename group column to group
    names(x)[names(x) == group] <- "group"

    # define magrittr pipe
    `%>%` <- magrittr::`%>%`

    # create the long file for the ICC routine/function default behavior is to sort alphabetically the else
    # routine keeps original variable order
    if (alpha.order) {
        long.dat <- x %>% tidyr::gather(type, score, -group)
    } else {
        long.dat <- x %>% tidyr::gather(type, score, -group) %>%
          dplyr::mutate(type = factor(type, levels = names(subset(x,
            select = -group))))  # add this line to convert the key to a factor
    }

    # New tidyverse version of mlm.iccs
    aov_model <- function(df) {
        lmr.model <- lmerTest::lmer(score ~ 1 + (1 | group), data = df)
    }

    aov_test <- function(df) {
        lmr.model <- lmerTest::lmer(score ~ 1 + (1 | group), data = df)
        ll.test <- lmerTest::ranova(lmr.model)
    }

    # get the model estimates
    models <- long.dat %>%
      tidyr::nest(-type) %>%
      dplyr::mutate(aov_obj = purrr::map(data, aov_model),
                    summaries = purrr::map(aov_obj,
                    broom.mixed::tidy)) %>%
      tidyr::unnest(summaries, .drop = T) %>%
      dplyr::select(type, effect, estimate,term) %>%
      dplyr::filter(effect != "fixed") %>%
      dplyr::mutate(variance = estimate^2) %>%
      dplyr::select(-estimate,-effect) %>%
      tidyr::spread(term, variance) %>%
      dplyr::rename(group.var = `sd__(Intercept)`,
                    residual = sd__Observation) %>%
      dplyr::mutate(ICC = group.var/(group.var + residual)) %>%
      dplyr::mutate(ICC = DescTools::Format(ICC,digits = 2, leading = "drop"))

    # get the ranova LRTs (and remove warnings from broom.mixed)
    options(warn = -1)
    tests <- long.dat %>%
      tidyr::nest(-type) %>%
      dplyr::mutate(test_obj = purrr::map(data, aov_test),
                    test_summaries = purrr::map(test_obj,
                    broom.mixed::tidy)) %>%
      tidyr::unnest(test_summaries, .drop = T) %>%
      dplyr::filter(!is.na(LRT))
    options(warn = 0)

    mlm.iccs <- long.dat %>%
      dplyr::group_by(type) %>%
      na.omit %>%
      dplyr::summarise(mean = mean(score, na.rm = T),
        sd = sd(score, na.rm = T), n = n()/dplyr::n_distinct(group)) %>%
      # dataframe needed because descTools chokes without it
        as.data.frame() %>%
        dplyr::left_join(., models, by = "type") %>%
        dplyr::left_join(., tests[c("type", "p.value")], by = "type") %>%
        dplyr::mutate(ICC2 = group.var/(group.var + residual/n)) %>%
        dplyr::select(type, mean, sd, ICC, p.value, ICC2) %>%
        dplyr::mutate(mean = DescTools::Format(mean, digits = 2,
                                               leading = "drop")) %>%
        dplyr::mutate(sd = DescTools::Format(sd, digits = 2,
                                             leading = "drop")) %>%
        dplyr::mutate(ICC2 = DescTools::Format(ICC2,digits = 2,
                                               leading = "drop")) %>%
        dplyr::mutate(icc.stars = ifelse(p.value < 0.01, paste0(ICC, "**"),
                                         ifelse(p.value < 0.05,
                                                paste0(ICC,"*"),ICC))) %>%
        dplyr::select(type, mean, sd, icc.stars, ICC2)

    mlm.iccs
    # Adapted from the corstars function https://github.com/Cogitos/statxp/blob/master/R/corstars.R

    # select the variables to be correlated
    cor.vars <- dplyr::select(x, -c(group))

    # sort columns by alpha order if indicated
    if (alpha.order) {
        cor.vars <- cor.vars %>% dplyr::select(sort(names(.)))
    }


    # Compute correlation matrix
    cor.vars <- as.matrix(cor.vars)
    correlation_matrix <- Hmisc::rcorr(cor.vars, type = "pearson")
    R <- correlation_matrix$r  # Matrix of correlation coeficients
    # print(R)
    p <- correlation_matrix$P  # Matrix of p-value

    # gets rid of the four stars bullshit
    mystars <- ifelse(p < 0.01, "**  ", ifelse(p < 0.05, "*   ", "    "))

    # get group mean centered correlations
    if (gmc) {

        cor.vars.c <- x

        # rename group variable
        cor.vars.c <- cor.vars.c %>% dplyr::rename(group.c = group)

        # sort columns by alpha order if indicated
        if (alpha.order) {
            cor.vars.c <- cor.vars.c %>% dplyr::select(sort(names(.)))
        }


        # group mean center the variables
        cor.vars.c[colnames(cor.vars.c)] <- lapply(cor.vars.c[colnames(cor.vars.c)],
                                                   function(y) y - ave(y,
                                                  cor.vars.c$group.c, FUN = mean))

        cor.vars.c <- dplyr::select(cor.vars.c, -c(group.c))

        correlation_matrix.c <- Hmisc::rcorr(as.matrix(cor.vars.c,
                                                       type = "pearson"))

        R.c <- correlation_matrix.c$r  # Matrix of correlation coeficients
        p.c <- correlation_matrix.c$P  # Matrix of p-value

        # gets rid of the four stars bullshit
        mystars.c <- ifelse(p.c < 0.01, "**  ", ifelse(p.c < 0.05, "*   ", "    "))

        mystars.c[lower.tri(mystars.c, diag = TRUE)] <- ""
        mystars[upper.tri(mystars, diag = TRUE)] <- ""

        R.c[lower.tri(R.c, diag = T)] <- 0
        R[upper.tri(R, diag = T)] <- 0

        R <- R + R.c

        R <- DescTools::Format(R, digits = 2, leading = "drop")

        all.stars <- matrix(paste0(mystars, mystars.c), ncol = ncol(R.c))

        ## build a new matrix that includes the correlations with their apropriate stars
        Rnew <- matrix(paste0(R, all.stars), ncol = ncol(R.c))
        diag(Rnew) <- paste0(diag(Rnew), "")
        rownames(Rnew) <- colnames(cor.vars.c)

        footer <- "<i>Note</i>: *<i>p</i> < .05. **<i>p</i> < .01. Correlations on the lower diagonal are at the individual level of analysis. Correlations on the upper diagonal are group-mean centered."

    } else {

        # easy way to get 2 decimals and drop leading 0
        R <- DescTools::Format(R, digits = 2, leading = "drop")

        ## build a new matrix that includes the correlations with their apropriate stars
        Rnew <- matrix(paste(R, mystars, sep = ""), ncol = ncol(cor.vars))
        diag(Rnew) <- paste(diag(R), "", sep = "")
        rownames(Rnew) <- colnames(cor.vars)

        ## remove upper triangle of correlation matrix
        Rnew <- as.matrix(Rnew)
        Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
        Rnew <- as.data.frame(Rnew, stringsAsFactors = FALSE)

        # Get column numbers as names
        col.nums <- rep(1:length(Rnew), 1)
        colnames(Rnew) <- as.character(col.nums)

        # row names preceded by number and period per apa
        row.nums <- rep(1:length(Rnew), 1)
        rownames(Rnew) <- paste(row.nums, ". ", rownames(Rnew), sep = "")

        footer <- "<i>Note</i>: *<i>p</i> < .05. **<i>p</i> < .01. Correlations are at the individual level of analysis."

    }

    # kill the diagonal
    diag(Rnew) <- "--"

    # get names for table
    tablenames <- as.character(colnames(cor.vars))
    # capitalize first letter in each variable name
    tablenames <- paste0(toupper(substr(tablenames, 1, 1)), substr(tablenames,
                                                            2, nchar(tablenames)))

    # bind the correlation tables
    cbind(mlm.iccs[-1], Rnew)

    if(result=="html") {
    # htmlTable
    htmlTable::htmlTable(cbind(mlm.iccs[-1], Rnew),
                         header = c("Mean", "SD", "ICC(1)",
                                    "ICC(2)", rep(1:ncol(Rnew))),
        rnames = paste(1:nrow(Rnew), ". ", tablenames, sep = ""),
        css.cell = "padding-left: 1em; padding-right: 1em;",
        caption = paste0("<b>", title, "</b>"), tfoot = footer)
} else {
  cbind(mlm.iccs[-1], Rnew)
}

    #Proof of concept


    # end icc.corrs function
}


#' Corstars
#'
#' This function creates the ICC matrix
#' @param x Data object.
#' @param method Correlation method. Default is pearson
#' @param removeTriangle Default is upper (per APA).
#' @param alpha.order Alphabetize variables.  Default is FALSE.
#' @param result Output options.  Default is html table to viewer.  Option "text" returns output to console only.
#' @param title Table caption.
#' @return A correlation table
#' @export


corstars <-function(x, method="pearson", removeTriangle=c("upper", "lower"), alpha.order = F,
                    result="html", sumstats=T, title="Correlation Table"){

  #define magrittr pipe
  `%>%` <- magrittr::`%>%`

  #sort columns by alpha order if requested
  #default is not
  if (alpha.order) {
    x <- x %>%
      dplyr::select(sort(names(.)))
  }

  #create duplicate of data frame for summary stats
  if(sumstats) {
    tempdf <- x
  }

  #Compute correlation matrix
  x <- as.matrix(x)
  correlation_matrix<-Hmisc::rcorr(x, type=method[1])
  R <- correlation_matrix$r # Matrix of correlation coeficients
  p <- correlation_matrix$P # Matrix of p-value

  ## Define notions for significance levels; spacing is important.
  #mystars <- ifelse(p < .0001, "****", ifelse(p < .001, "*** ", ifelse(p < .01, "**  ", ifelse(p < .05, "*   ", "    "))))

  #gets rid of the four stars bullshit
  mystars <- ifelse(p < .01, "**  ", ifelse(p < .05, "*   ", "    "))


  ## trunctuate the correlation matrix to two decimal
  #R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1]

  R <- DescTools::Format(R, digits=2,leading="drop", na.form="--", sci = NA)
  #print(R)

  ## build a new matrix that includes the correlations with their apropriate stars
  Rnew <- matrix(paste0(R, mystars), ncol=ncol(x))
  diag(Rnew) <- paste0(diag(R), "")
  rownames(Rnew) <- colnames(x)


  #colnames(Rnew) <- paste(colnames(x), "", sep="")

  ## remove upper triangle of correlation matrix
  if(removeTriangle[1]=="upper"){
    Rnew <- as.matrix(Rnew)
    Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
    Rnew <- as.data.frame(Rnew, stringsAsFactors = F)
  }

  ## remove lower triangle of correlation matrix
  else if(removeTriangle[1]=="lower"){
    Rnew <- as.matrix(Rnew)
    Rnew[lower.tri(Rnew, diag = TRUE)] <- ""
    Rnew <- as.data.frame(Rnew, stringsAsFactors = F)
  }


  #Get column numbers as names
  col.nums <- rep(1:length(Rnew),1)
  colnames(Rnew) <- as.character(col.nums)
  diag(Rnew) <- "--"

  if(sumstats) {
    #get just the mean and SD
    tempstats <- as.data.frame(psych::describe(tempdf))[3:4]
    tempstats <- DescTools::Format(tempstats, digits=2, na.form="")

    #insert descriptive stats at front of table
    Rnew <- cbind(tempstats, Rnew)
    #process htmltable info
    table.headers <- c("Mean","SD",rep(1:ncol(R)))
    tempcol <- rep("", nrow(Rnew))
    cell.form <- cbind(tempcol,tempcol,R)
  } else {
    table.headers <- rep(1:ncol(R))
  }


  if (result[1]=="text")  {
    rownames(Rnew) <- paste0(toupper(substr(rownames(Rnew), 1, 1)),
                             substr(rownames(Rnew), 2, nchar(rownames(Rnew))))
    row.nums <- rep(1:length(rownames(Rnew)),1)

    rownames(Rnew) <- paste(row.nums,". ", rownames(Rnew), sep = "")

    return(Rnew)

  } else if (result[1]=="html") {

    #get names for table
    tablenames <- as.character(colnames(x))
    #capitalize first letter in each variable name
    tablenames <- paste0(toupper(substr(tablenames, 1, 1)), substr(tablenames, 2, nchar(tablenames)))

    htmlTable::htmlTable(Rnew, header=table.headers,
                         rnames = paste0(1:nrow(Rnew), ". ", tablenames),
                         css.cell = "padding-left: .5em; padding-right: .2em;",
                         caption = title,
                         rowlabel="Variables",
                         tfoot="<i>Note</i>: *<i>p</i> < .05. **<i>p</i> < .01.")
  } else {
    xtable::xtable(Rnew, type="latex")
  }
  #ends function
}

#' Latent Group Model
#'
#' This function creates the LGM matrix
#' @param x Data object.
#' @param group Nesting variable.
#' @param title Table capation.
#' @param stars Provide significance stars.  Default is TRUE.
#' @param sumstats Provide summary stats.  Default is TRUE.
#' @param sumtable Provide sem summary table. Default is FALSE.
#' @param alpha.order Sort variables alphabetically.  Default is FALSE.
#' @return A correlation table
#' @export

lgm <-function(x, group, title="LGM", stars=TRUE, sumstats=TRUE, sumtable = FALSE, alpha.order = FALSE) {
  options(scipen=999)

  #define magrittr pipe
  `%>%` <- magrittr::`%>%`

  #ungroup if original dataset is grouped
  if(dplyr::is.grouped_df(x)){

    x <- x %>%
      dplyr::ungroup()
  }

  #remove group variable from df
  tempdf <- dplyr::select(x, -c(group))

  #sort columns by alpha order if indicated
  if(alpha.order) {
  tempdf <- tempdf %>%
    dplyr::select(sort(names(.)))
  }

  #get number of variables
  varnum <- ncol(tempdf)

  #get variable names in a string
  varnames <- paste( unlist(colnames(tempdf)), collapse=' ')

  #message(print(str(varnames)))

  #get names for table
  tablenames <- as.character(colnames(tempdf))
  #capitalize first letter in each variable name
  tablenames <- paste0(toupper(substr(tablenames, 1, 1)), substr(tablenames, 2, nchar(tablenames)))

  #initialize loop variables
  dep <- NULL
  inds <- NULL
  equation <- NULL

  #set counter
  y=1

  while (y < varnum) {
    #message(print(varnames))


    #set up left side of equation
    dep <- sub(" +","~~", varnames)

    #set up right side of equation
    inds <- gsub(" +","\\+", dep)

    #paste together equation through loop
    equation <- paste(equation,inds,"\n",sep="")

    #prep varnames for next iteration through loop
    #this removes first word/variables
    varnames <- stringr::word(varnames, 2, -1)


    #reset dep and increment counter
    dep <- NULL
    y <- y+1
  }

  #put lavaan text as needed
  model.eq <- paste0("level: 1\n", equation, "level: 2\n", equation, sep="")

  suppressMessages(library(lavaan))

  model <- model.eq

  model.out <- sem(model, data = x, cluster = group,
                   std.lv = F, verbose = FALSE)

  if(sumtable) {
    summary(model.out, standardized=TRUE)
  }

  #now build and output LGM matrix
  #correlations printed but significance for covariances
  ind <- cov2cor(as.matrix(as.data.frame(inspect(model.out,"cov.ov")[[1]])))
  grp <- cov2cor(as.matrix(as.data.frame(inspect(model.out,"cov.ov")[[2]])))

  all.corrs.sem <- ind*0
  ind[upper.tri(ind, diag = TRUE)] <- 0
  grp[lower.tri(grp, diag = TRUE)] <- 0
  all.corrs.sem <- ind + grp
  diag(all.corrs.sem) <- diag(as.matrix(as.data.frame(inspect(model.out,"cov.ov")[[2]])))/
    (diag(as.matrix(as.data.frame(inspect(model.out,"cov.ov")[[2]])))+
       diag(as.matrix(as.data.frame(inspect(model.out,"cov.ov")[[1]]))))

  #create the cell formatting for htmltable
  cell.form <- all.corrs.sem*0
  diag(cell.form) <- rep("font-weight: bold",nrow(cell.form))
  cell.form[upper.tri(cell.form, diag=FALSE)] <- ""
  cell.form[lower.tri(cell.form, diag=FALSE)] <- "font-style: italic;"

  #removes leading zeros and rounds to two decimal places
  all.corrs.sem <- DescTools::Format(all.corrs.sem, digits=2,leading="drop", na.form="--")

  #include stars by default
  if(stars) {
    z.w.groups <- lavInspect(model.out, "est")$within$theta/lavInspect(model.out, "se")$within$theta
    z.b.groups <- as.matrix(as.data.frame(inspect(model.out, "est")[[2]]["theta"]))/as.matrix(as.data.frame(inspect(model.out, "se")[[2]]["theta"]))

    #get table of probabilities
    #diag is for the variances
    prob.inds <- pnorm(z.w.groups, lower.tail = FALSE)*2
    prob.groups <- pnorm(z.b.groups, lower.tail = FALSE)*2

    #gets rid of the four stars bullshit
    mystars.ind <- ifelse(prob.inds < .01, "**  ", ifelse(prob.inds < .05, "*   ", "    "))
    mystars.grp <- ifelse(prob.groups < .01, "**  ", ifelse(prob.groups < .05, "*   ", "    "))

    mystars.grp[lower.tri(mystars.grp, diag=FALSE)] <- ""
    mystars.ind[upper.tri(mystars.ind, diag=TRUE)] <- ""

    all.stars <- matrix(paste0(mystars.ind, mystars.grp), ncol=ncol(prob.inds))

    matrix.out <- matrix(paste0(all.corrs.sem, all.stars), ncol=ncol(all.corrs.sem))
    table.footer <- "<i>Note</i>: *<i>p</i> < .05. **<i>p</i> < .01. Intraclass correlations are on the diagonal (in bold). Individual-level correlations are on the lower diagonal (in italics). Group-level correlations are on the upper diagonal."
  } else {
    matrix.out <- all.corrs.sem
    table.footer <- "Intraclass correlations are on the diagonal (in bold). Individual-level correlations are on the lower diagonal (in italics). Group-level correlations are on the upper diagonal."
  }

  if(sumstats) {
    #get just the mean and SD
    tempstats <- as.data.frame(psych::describe(tempdf))[3:4]
    tempstats <- DescTools::Format(tempstats, digits=2, na.form="")

    #insert descriptive stats at front of table
    matrix.out <- cbind(tempstats, matrix.out)
    #process htmltable info
    table.headers <- c("Mean","SD",rep(1:ncol(all.corrs.sem)))
    tempcol <- rep("", nrow(all.corrs.sem))
    cell.form <- cbind(tempcol,tempcol, cell.form)
  } else {
    table.headers <- rep(1:ncol(all.corrs.sem))
  }


  library(htmlTable)
  htmlTable(matrix.out,
            header=table.headers,
            rnames = paste(1:nrow(all.corrs.sem), ". ", tablenames, sep = ""),
            #css.cell = "padding-left: 1em; padding-right: 1em;",
            css.cell = cell.form,
            caption = paste0("<b>",title,"</b>"),
            tfoot = table.footer)
  #ends the lgm function
}


