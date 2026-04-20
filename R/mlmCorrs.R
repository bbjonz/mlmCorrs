#' Estimate ICCs, correlations, and descriptive Statistics
#'
#' This function creates the ICC matrix
#' @param x Data object.
#' @param group Nesting variable.
#' @param title Table caption.
#' @param gmc Provide group-mean center correlations.  Default is FALSE.
#' @param stars Number of significance stars.  Default is 2, max is 4 (p < .0001)
#' @param alpha.order Sort variables alphabetically.  Default is False
#' @param result Output options.  Default is kable to viewer.  Option "text" returns output to console only.
#' @param icc2 Report ICC2 (Bliese).  Default is not
#' @return A correlation table with sample statistics and ICC estimates
#' @export

# Descriptive Stats with ICCs and Correlations ####
icc.corrs <- function(x, group, title = "Descriptive Stats",
                      gmc = FALSE, stars = 2,
                      alpha.order = FALSE, result = "html",
                      icc2 = FALSE) {

  list.of.packages <- c("tidyverse","psych","lme4","Hmisc","DescTools","knitr","kableExtra","gt")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)

  options(scipen=999)

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
        long.dat <- x %>% tidyr::pivot_longer(!group, names_to="type", values_to="score")
    } else {
        long.dat <- x %>% tidyr::pivot_longer(!group, names_to="type", values_to="score") %>%
          dplyr::mutate(type = factor(type, levels = names(subset(x,
            select = -group))))  # add this line to convert the key to a factor
    }

    # New tidyverse version of mlm.iccs
    # Get model estimates
    aov_model <- function(df) {
        lmr.model <- lmerTest::lmer(score ~ 1 + (1 | group), data = df)
    }

    #Get random effect significance test
    aov_test <- function(df) {
        lmr.model <- lmerTest::lmer(score ~ 1 + (1 | group), data = df)
        ll.test <- lmerTest::ranova(lmr.model)
    }

    #get sample size for both levels
    aov_model.ss <- function(df) {
      lmr.model <- lmerTest::lmer(score ~ 1 + (1 | group), data = df)
      size = c(NA, lmr.model@pp$Ut@Dim[1], lmr.model@pp$Ut@Dim[2])
      return(data.frame(size))
    }

    # get the model estimates
    models <- long.dat %>%
      tidyr::nest(-type) %>%
      dplyr::mutate(aov_obj = purrr::map(data, aov_model),
                    sampsize = purrr::map(data, aov_model.ss),
                    summaries = purrr::map(aov_obj,
                    broom.mixed::tidy)) %>%
      tidyr::unnest(c(summaries, sampsize)) %>%
      dplyr::select(type, effect, estimate, group, size) %>%
      dplyr::filter(effect != "fixed") %>%
      dplyr::mutate(variance = estimate^2) %>%
      dplyr::select(-estimate,-effect) %>%
      tidyr::pivot_wider(id_cols=type, values_from = c(variance, size),
                         names_from = group) %>%
      dplyr::rename(group.var = variance_group,
                    residual = variance_Residual,
                    "N Groups" = size_group,
                    "Total N" = size_Residual) %>%
      dplyr::mutate(ICC = group.var/(group.var + residual)) %>%
      dplyr::mutate(ICC = sub("^(-?)0.", "\\1.", sprintf("%.2f", ICC)))

    # get the ranova LRTs (and remove warnings from broom.mixed)
    options(warn = -1)
    tests <- long.dat %>%
      tidyr::nest(-type) %>%
      dplyr::mutate(test_obj = purrr::map(data, aov_test),
                    test_summaries = purrr::map(test_obj,
                    broom.mixed::tidy)) %>%
      tidyr::unnest(test_summaries) %>%
      dplyr::filter(!is.na(LRT))
    options(warn = 0)

    mlm.iccs <- long.dat %>%
      dplyr::group_by(type) %>%
      na.omit %>%
      dplyr::summarise(mean = mean(score, na.rm = T),
        sd = sd(score, na.rm = T), n = n()/dplyr::n_distinct(group)) %>%
        as.data.frame() %>%
        dplyr::left_join(., models, by = "type") %>%
        dplyr::left_join(., tests[c("type", "p.value")], by = "type") %>%
        dplyr::mutate(ICC2 = group.var/(group.var + residual/n)) %>%
        dplyr::select(type, "N Groups", "Total N", mean, sd, ICC, p.value, ICC2) %>%
        dplyr::mutate(mean = sub("^(-?)0.", "\\1.", sprintf("%.2f", mean))) %>%
        dplyr::mutate(sd = sub("^(-?)0.", "\\1.", sprintf("%.2f", sd))) %>%
        dplyr::mutate(ICC2 = sub("^(-?)0.", "\\1.", sprintf("%.2f", ICC2))) %>%
        dplyr::mutate(ICC1 = ifelse(p.value < 0.01, paste0(ICC, "**"),
                                         ifelse(p.value < 0.05,
                                                paste0(ICC,"*"),ICC))) %>%
        dplyr::select(type, "N Groups", "Total N", mean, sd, ICC1, ICC2)

    #print(mlm.iccs)
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

    if(stars == 2) {
      mystars <- ifelse(p < .01, "**  ", ifelse(p < .05, "*   ", "    "))
      #footer for table
      footer <- "*<i>p</i> < .05. **<i>p</i> < .01."
    } else if(stars == 3) {
      #footer for table
      mystars <- ifelse(p < .001, "*** ", ifelse(p < .01, "**  ", ifelse(p < .05, "*   ", "    ")))
      footer <- "*<i>p</i> < .05. **<i>p</i> < .01. ***<i>p<i/> < .001. "
    } else if(stars == 4) {
      mystars <- ifelse(p < .0001, "****", ifelse(p < .001, "*** ", ifelse(p < .01, "**  ", ifelse(p < .05, "*   ", "    "))))
      footer <- "*<i>p</i> < .05. **<i>p</i> < .01. ***<i>p</i> < .001. ****<i>p</i> < .0001."

    } else {
      stop("You requested more than 4 significance stars.  Please provide a valid number between 2 and 4")
    }

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

        R.c <- correlation_matrix.c$r  # Matrix of correlation coefficients
        p.c <- correlation_matrix.c$P  # Matrix of p-value

        if(stars == 2) {
          mystars.c <- ifelse(p.c < .01,
                              "**  ",
                              ifelse(p.c < .05, "*   ", "    "))
          #footer for table
          footer <- "*<i>p</i> < .05. **<i>p</i> < .01."
        } else if(stars == 3) {
          #footer for table
          mystars.c <- ifelse(p.c < .001, "*** ",
                              ifelse(p.c < .01,"**  ",
                                     ifelse(p.c < .05, "*   ", "    ")))
          footer <- "*<i>p</i> < .05. **<i>p</i> < .01. ***<i>p<i/> < .001. "
        } else if(stars == 4) {
          mystars.c <- ifelse(p.c < .0001, "****",
                              ifelse(p.c < .001, "*** ",
                                     ifelse(p.c < .01, "**  ",
                                            ifelse(p.c < .05, "*   ", "    "))))
          footer <- "*<i>p</i> < .05. **<i>p</i> < .01. ***<i>p</i> < .001. ****<i>p</i> < .0001."

        } else {
          stop("You requested more than 4 significance stars.  Please provide a valid number between 2 and 4")
        }

        mystars.c[lower.tri(mystars.c, diag = TRUE)] <- ""
        mystars[upper.tri(mystars, diag = TRUE)] <- ""

        R.c[lower.tri(R.c, diag = T)] <- 0
        R[upper.tri(R, diag = T)] <- 0

        R <- R + R.c

        #R <- DescTools::Format(R, digits = 2, leading = "drop")
        R <- matrix(sub("^(-?)0.", "\\1.", sprintf("%.2f", R)), nrow = nrow(R))

        all.stars <- matrix(paste0(mystars, mystars.c), ncol = ncol(R.c))

        ## build a new matrix that includes the correlations with their appropriate stars
        Rnew <- matrix(paste0(R, all.stars), ncol = ncol(R.c))
        diag(Rnew) <- paste0(diag(Rnew), "")
        rownames(Rnew) <- colnames(cor.vars.c)

        footer <- "* <i>p</i> < .05. **<i>p</i> < .01. Correlations on the lower diagonal are at the individual level of analysis. Correlations on the upper diagonal are group-mean centered."

    } else {

        # easy way to get 2 decimals and drop leading 0
        R <- matrix(sub("^(-?)0.", "\\1.", sprintf("%.2f", R)), nrow = nrow(R))

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

        footer <- "*<i>p</i> < .05. **<i>p</i> < .01. Correlations are at the individual level of analysis."

    }

    # kill the diagonal
    diag(Rnew) <- "--"

    # get names for table
    tablenames <- as.character(colnames(cor.vars))
    # capitalize first letter in each variable name
    tablenames <- paste0(toupper(substr(tablenames, 1, 1)), substr(tablenames,
                                                            2, nchar(tablenames)))

    # bind the correlation tables and provide row/col names
    # ICC2 not reported by default
    cbind(mlm.iccs[-1], Rnew)
    tablePrint <- cbind(mlm.iccs[-1], Rnew)
    row.names(tablePrint) <- paste0(1:nrow(Rnew), ". ", tablenames)
    colnames(tablePrint) <- c("N Groups", "Total N", "Mean", "SD", "ICC(1)",
                              "ICC(2)", rep(1:ncol(Rnew)))

    if (!icc2) {
      tablePrint %>%
        dplyr::select(-"ICC(2)") %>%
        dplyr::rename(ICC = "ICC(1)") -> tablePrint
    }

    if(result=="html") {
      tablePrint %>%
        tibble::rownames_to_column(.,"Variable") %>%
        gt::gt(caption = title) %>%
        gt::tab_options(table.border.bottom.width = "0px",
                        table.border.top.width = "0px",
                        heading.align = "left") %>%
        #gt::tab_header(title = title) %>%
        gt::tab_source_note(
          source_note = gt::html(c("<i>Note</i>. ", footer))
        )

    } else if (result[1]=="text") {
    return(tablePrint)
} else {
  cbind(mlm.iccs[-1], Rnew)
}

    # ends the icc.corrs function
}

# APA Correlation Table ####
#' Corstars
#'
#' This function creates a correlation matrix with descriptive statistics.
#' @param x Data object.
#' @param group Optional grouping variable (unquoted column name) for stratified tables.
#' @param method Correlation method. Default is pearson
#' @param removeTriangle Default is upper (per APA).
#' @param alpha.order Alphabetize variables. Default is FALSE.
#' @param stars Number of significance stars. Default is 2, max is 4 (p < .0001)
#' @param result Output options. Default is "html". Option "text" returns a data frame.
#' @param sumstats Include mean, SD, and N. Default is TRUE.
#' @param title Table caption.
#' @return A correlation table (gt object, data frame, or xtable)
#' @importFrom rlang enquo quo_is_null quo_name
#' @export

corstars <- function(x, group = NULL, method = "pearson",
                     removeTriangle = c("upper", "lower"),
                     alpha.order = FALSE, stars = 2, result = "html",
                     sumstats = TRUE, title = "Correlation Table") {

  list.of.packages <- c("tidyverse", "psych", "Hmisc", "DescTools",
                        "knitr", "kableExtra", "gt", "rlang")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if (length(new.packages)) install.packages(new.packages)

  options(scipen = 999)
  `%>%` <- magrittr::`%>%`

  # ── capture group argument safely for package use ─────────────────────────
  group_quo <- rlang::enquo(group)
  is_grouped <- !rlang::quo_is_null(group_quo)
  group_var  <- if (is_grouped) rlang::quo_name(group_quo) else NULL

  # ── footer ────────────────────────────────────────────────────────────────
  footer <- switch(as.character(stars),
                   "2" = "*<i>p</i> < .05. **<i>p</i> < .01.",
                   "3" = "*<i>p</i> < .05. **<i>p</i> < .01. ***<i>p</i> < .001.",
                   "4" = "*<i>p</i> < .05. **<i>p</i> < .01. ***<i>p</i> < .001. ****<i>p</i> < .0001.",
                   stop("Please provide a valid number of stars between 2 and 4.")
  )

  # ── helper: capitalize first letter only ─────────────────────────────────
  cap_first <- function(s) {
    paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
  }

  # ── shared gt formatting ──────────────────────────────────────────────────
  format_gt <- function(gt_obj) {
    gt_obj %>%
      gt::tab_options(
        table.border.bottom.width = "0px",
        table.border.top.width    = "0px",
        heading.align             = "left",
        row_group.font.weight     = "bold"
      ) %>%
      gt::tab_header(title = title) %>%
      gt::cols_align(align = "center", columns = everything()) %>%
      gt::cols_align(align = "left",   columns = "Variable") %>%
      gt::tab_source_note(
        source_note = gt::html(c("<i>Note</i>. ", footer))
      )
  }

  # ── inner workhorse: builds one formatted data frame for a data slice ─────
  build_table <- function(df, group_label = NULL) {

    if (alpha.order) df <- df %>% dplyr::select(sort(names(.)))

    tempdf <- df

    # Correlation matrix
    x_mat <- as.matrix(df)
    cm    <- Hmisc::rcorr(x_mat, type = method[1])
    R     <- cm$r
    p     <- cm$P
    ntemp <- cm$n

    # Significance stars
    mystars <- switch(as.character(stars),
                      "2" = ifelse(p < .01,  "**  ", ifelse(p < .05, "*   ", "    ")),
                      "3" = ifelse(p < .001, "*** ", ifelse(p < .01, "**  ",
                                                            ifelse(p < .05, "*   ", "    "))),
                      "4" = ifelse(p < .0001,"****", ifelse(p < .001, "*** ",
                                                            ifelse(p < .01, "**  ", ifelse(p < .05, "*   ", "    "))))
    )

    # Format r values (drop leading zero)
    R_fmt <- matrix(sub("^(-?)0.", "\\1.", sprintf("%.2f", R)), nrow = nrow(R))

    # Paste r + stars
    Rnew <- matrix(paste0(R_fmt, mystars), ncol = ncol(x_mat))
    diag(Rnew) <- paste0(diag(R_fmt), "")
    rownames(Rnew) <- colnames(x_mat)

    # Remove triangle
    if (removeTriangle[1] == "upper") {
      Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
    } else {
      Rnew[lower.tri(Rnew, diag = TRUE)] <- ""
    }
    Rnew <- as.data.frame(Rnew, stringsAsFactors = FALSE)

    # Column numbers as names; diagonal = "--"
    colnames(Rnew) <- as.character(seq_len(ncol(Rnew)))
    diag(Rnew) <- "--"

    # Summary stats
    if (sumstats) {
      tempstats <- data.frame(
        Mean = colMeans(tempdf, na.rm = TRUE),
        SD   = apply(tempdf, 2, sd, na.rm = TRUE),
        N    = diag(ntemp)
      )
      tempstats <- as.data.frame(lapply(tempstats, sprintf, fmt = "%.2f"))
      Rnew <- cbind(tempstats, Rnew)
      names(Rnew) <- c("Mean", "SD", "N", seq_len(ncol(R)))
    }

    # Numbered row names
    rownames(Rnew) <- paste0(
      seq_len(nrow(Rnew)), ". ",
      cap_first(rownames(Rnew))
    )

    attr(Rnew, "group_label") <- group_label
    Rnew
  }

  # ── build tables ──────────────────────────────────────────────────────────
  if (!is_grouped) {
    # ── No grouping: single table ───────────────────────────────────────────
    numeric_cols <- x %>% dplyr::select(where(is.numeric))
    tables <- list(build_table(numeric_cols, group_label = NULL))

  } else {
    # ── Grouped: one table per level ────────────────────────────────────────
    if (!group_var %in% names(x))
      stop(paste0("Grouping variable '", group_var, "' not found in data."))

    grp_vec  <- x[[group_var]]
    grp_lvls <- sort(unique(grp_vec[!is.na(grp_vec)]))

    tables <- lapply(grp_lvls, function(lvl) {
      sub_df <- x %>%
        dplyr::filter(.data[[group_var]] == lvl) %>%
        dplyr::select(where(is.numeric))
      build_table(sub_df, group_label = cap_first(as.character(lvl)))
    })
  }

  # ── return results ────────────────────────────────────────────────────────
  if (result[1] == "text") {
    if (length(tables) == 1) return(tables[[1]])

    do.call(rbind, lapply(seq_along(tables), function(i) {
      tbl <- tables[[i]]
      lbl <- attr(tbl, "group_label")
      sep <- as.data.frame(
        matrix("", nrow = 1, ncol = ncol(tbl)),
        stringsAsFactors = FALSE
      )
      names(sep)    <- names(tbl)
      rownames(sep) <- paste0("── ", lbl, " ──")
      rbind(sep, tbl)
    }))

  } else if (result[1] == "html") {

    if (!is_grouped) {
      # ── Single table: no groupname_col, no spurious header ──────────────
      tables[[1]] %>%
        tibble::rownames_to_column("Variable") %>%
        gt::gt() %>%
        format_gt()

    } else {
      # ── Grouped: attach label column and use groupname_col ───────────────
      all_rows <- lapply(tables, function(tbl) {
        df <- tibble::rownames_to_column(tbl, "Variable")
        df$`.group` <- attr(tbl, "group_label")
        df
      })

      do.call(rbind, all_rows) %>%
        gt::gt(groupname_col = ".group") %>%
        format_gt()
    }

  } else {
    # LaTeX
    combined_list <- lapply(tables, function(tbl) {
      lbl <- attr(tbl, "group_label")
      if (!is.null(lbl)) cat("\n%% Group:", lbl, "\n")
      xtable::xtable(tbl)
    })
    invisible(combined_list)
  }
}

# Latent Group Model ####
#' This function creates the LGM matrix
#' @param x Data object.
#' @param group Nesting variable.
#' @param title Table caption.
#' @param printstars Provide significance stars.  Default is TRUE.
#' @param result Output options.  Default is html table to viewer.  Option "text" returns just that.
#' @param sumstats Provide summary stats.  Default is TRUE.
#' @param sumtable Provide sem summary table. Default is FALSE.
#' @param stars Number of significance stars.  Default is 2, max is 4 (p < .0001)
#' @param alpha.order Sort variables alphabetically.  Default is FALSE.
#' @return A correlation table
#' @export

lgm <-function(x, group, title="LGM", printstars=TRUE, result = "html",
               sumstats=TRUE, sumtable = FALSE, stars = 2, alpha.order = FALSE) {

  options(scipen=999)

  list.of.packages <- c("tidyverse","lavaan","psych","Hmisc","DescTools","knitr","kableExtra")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)

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
  ind <- cov2cor(as.matrix(as.data.frame(lavInspect(model.out,"cov.ov")[[1]])))
  grp <- cov2cor(as.matrix(as.data.frame(lavInspect(model.out,"cov.ov")[[2]])))

  all.corrs.sem <- ind*0
  ind[upper.tri(ind, diag = TRUE)] <- 0
  grp[lower.tri(grp, diag = TRUE)] <- 0
  all.corrs.sem <- ind + grp
  diag(all.corrs.sem) <- diag(as.matrix(as.data.frame(lavInspect(model.out,"cov.ov")[[2]])))/
    (diag(as.matrix(as.data.frame(lavInspect(model.out,"cov.ov")[[2]])))+
       diag(as.matrix(as.data.frame(lavInspect(model.out,"cov.ov")[[1]]))))

  #removes leading zeros and rounds to two decimal places
  all.corrs.sem <- matrix(sub("^(-?)0.", "\\1.", sprintf("%.2f", all.corrs.sem)), nrow = nrow(all.corrs.sem))

  rownames(all.corrs.sem) <- rownames(ind)
  colnames(all.corrs.sem) <- colnames(ind)

  #include stars by default
  if(printstars) {
    z.w.groups <- lavInspect(model.out, "est")$within$theta/lavInspect(model.out, "se")$within$theta
    z.b.groups <- as.matrix(as.data.frame(lavInspect(model.out, "est")[[2]]["theta"]))/as.matrix(as.data.frame(lavInspect(model.out, "se")[[2]]["theta"]))

    #get table of probabilities
    #diag is for the variances
    prob.inds <- pnorm(z.w.groups, lower.tail = FALSE)*2
    prob.groups <- pnorm(z.b.groups, lower.tail = FALSE)*2

    if(stars == 2) {
      mystars.ind <- ifelse(prob.inds < .01, "**  ",
                            ifelse(prob.inds < .05, "*   ", "    "))
      mystars.grp <- ifelse(prob.groups < .01, "**  ",
                            ifelse(prob.groups < .05, "*   ", "    "))
      #footer for table
      footer <- "*<i>p</i> < .05. **<i>p</i> < .01."
    } else if(stars == 3) {
      mystars.ind <- ifelse(prob.inds < .001, "*** ",
                            ifelse(prob.inds < .01, "**  ",
                                   ifelse(prob.inds < .05, "*   ", "    ")))
      mystars.grp <- ifelse(prob.groups < .001, "*** ",
                            ifelse(prob.groups < .01, "**  ",
                                   ifelse(prob.groups < .05, "*   ", "    ")))
      #footer for table
      footer <- "*<i>p</i> < .05. **<i>p</i> < .01. ***<i>p<i/> < .001. "
    } else if(stars == 4) {
      mystars.ind <- ifelse(prob.inds < .0001, "****",
                            ifelse(prob.inds < .001, "*** ",
                                   ifelse(prob.inds < .01, "**  ", ifelse(prob.inds < .05, "*   ", "    "))))
      mystars.grp <- ifelse(prob.groups < .0001, "****",
                            ifelse(prob.groups < .001, "*** ",
                                   ifelse(prob.groups < .01, "**  ", ifelse(prob.groups < .05, "*   ", "    "))))
      footer <- "*<i>p</i> < .05. **<i>p</i> < .01. ***<i>p</i> < .001. ****<i>p</i> < .0001."

    } else {
      stop("You requested more than 4 significance stars.  Please provide a valid number between 2 and 4")
    }

    mystars.grp[lower.tri(mystars.grp, diag=FALSE)] <- ""
    mystars.ind[upper.tri(mystars.ind, diag=TRUE)] <- ""

    all.stars <- matrix(paste0(mystars.ind, mystars.grp), ncol=ncol(prob.inds))

    matrix.out <- matrix(paste0(all.corrs.sem, all.stars), ncol=ncol(all.corrs.sem))

    #format diagonal with bold
    if (result=="html") {
    matrix.d <- diag(matrix.out)
    diag(matrix.out) <- matrix.d
    }

    table.footer <- "*<i>p</i> < .05. **<i>p</i> < .01. Intraclass correlations are on bold the diagonal. Individual-level correlations are on the lower diagonal. Group-level correlations are on the upper diagonal."
  } else {
    matrix.out <- all.corrs.sem
    table.footer <- "Intraclass correlations are on the diagonal. Individual-level correlations are on the lower diagonal. Group-level correlations are on the upper diagonal."
  }

  if(sumstats) {
    #get just the mean and SD
    tempstats <- data.frame(mean=colMeans(tempdf, na.rm = T), sd=apply(tempdf, 2, sd, na.rm = T))
    tempstats <- as.data.frame(lapply(tempstats, sprintf, fmt="%.2f"))

    #tempstats <- as.data.frame(psych::describe(tempdf))[3:4]
    #tempstats <- DescTools::Format(tempstats, digits=2, na.form="")

    #insert descriptive stats at front of table
    matrix.out <- cbind(tempstats, matrix.out)

    rownames(matrix.out) <- rownames(ind)
    #print(matrix.out)

    #process header info
    names(matrix.out) <- c("Mean","SD",rep(1:ncol(all.corrs.sem)))
}

    #process row names
    rownames(matrix.out) <- paste0(toupper(substr(rownames(matrix.out), 1, 1)),
                             substr(rownames(matrix.out), 2, nchar(rownames(matrix.out))))
    row.nums <- rep(1:length(rownames(matrix.out)),1)

    rownames(matrix.out) <- paste(row.nums,". ", rownames(matrix.out), sep = "")

    #note: gt package doesn't allow for bolding the diagonals yet.  Working on it.
    if (result=="html") {
      matrix.out %>%
        tibble::rownames_to_column(.,"Variable") %>%
        gt::gt () %>%
        gt::tab_options(table.border.bottom.width = "0px",
                        table.border.top.width = "0px",
                        heading.align = "left") %>%
        gt::tab_header(title = title) %>%
        gt::tab_source_note(
          source_note = gt::html(c("<i>Note</i>. ", footer))
        )


    } else {
      return(matrix.out)
    }
  #ends the lgm function
}


