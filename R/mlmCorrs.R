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


# APA Correlation Table Groups ####
#' Create an APA-style correlation table with optional group stacking
#'
#' @param x Data frame of numeric variables.
#' @param method Correlation method. Default is "pearson".
#' @param removeTriangle Which triangle to remove. Default is "upper" (APA).
#' @param alpha.order Alphabetize variables. Default is FALSE.
#' @param stars Number of significance stars. Default is 2, max is 4.
#' @param result Output format. Default is "html" (gt table). "text" returns
#'   the raw data frame.
#' @param sumstats Logical. If TRUE (default), include Mean, SD, and N.
#' @param title Character string for the table title.
#' @param group Optional character string naming a grouping column in x.
#'   If provided, a stacked gt table is produced with one section per group.
#' @param decimals Number of decimal places to display. Default is 2.
#'
#' @return A gt table object or data frame depending on result argument.

corstars <- function(x,
                          method = "pearson",
                          removeTriangle = c("upper", "lower"),
                          alpha.order = FALSE,
                          stars = 2,
                          result = "html",
                          sumstats = TRUE,
                          title = "Correlation Table",
                          group = NULL,
                          decimals = 2) {


  options(scipen = 999)

  # ── Significance star footer ───────────────────────────────────────────────
  if (stars == 2) {
    footer <- "*<i>p</i> < .05. **<i>p</i> < .01."
  } else if (stars == 3) {
    footer <- "*<i>p</i> < .05. **<i>p</i> < .01. ***<i>p</i> < .001."
  } else if (stars == 4) {
    footer <- "*<i>p</i> < .05. **<i>p</i> < .01. ***<i>p</i> < .001. ****<i>p</i> < .0001."
  } else {
    stop("Please provide a valid number of stars between 2 and 4.")
  }

  # ── Core helper: build one correlation block from a data frame ─────────────
  build_block <- function(df, grp_label = NULL) {

    if (alpha.order) df <- df |> dplyr::select(sort(names(df)))

    tempdf <- df

    # Correlation matrix
    x_mat <- as.matrix(df)
    corr   <- Hmisc::rcorr(x_mat, type = method[1])
    R      <- corr$r
    p      <- corr$P
    ntemp  <- corr$n

    # Significance stars
    mystars <- if (stars == 2) {
      ifelse(p < .01, "**  ", ifelse(p < .05, "*   ", "    "))
    } else if (stars == 3) {
      ifelse(p < .001, "*** ", ifelse(p < .01, "**  ",
                                      ifelse(p < .05, "*   ", "    ")))
    } else {
      ifelse(p < .0001, "****", ifelse(p < .001, "*** ",
                                       ifelse(p < .01, "**  ", ifelse(p < .05, "*   ", "    ") )))
    }

    # Format correlations — drop leading zero per APA and apply decimal formatting
    R_fmt <- R

    # Convert to character matrix and format each element
    R_char <- matrix("", nrow = nrow(R_fmt), ncol = ncol(R_fmt))

    for(i in 1:nrow(R_fmt)) {
      for(j in 1:ncol(R_fmt)) {
        if(!is.na(R_fmt[i,j])) {
          # Format to specified decimals with trailing zeros
          formatted_val <- sprintf(paste0("%.", decimals, "f"), R_fmt[i,j])
          # Remove leading zero if it exists and the value is between -1 and 1
          if(grepl("^-?0\\.", formatted_val) && abs(R_fmt[i,j]) < 1) {
            formatted_val <- sub("^(-?)0\\.", "\\1.", formatted_val)
          }
          R_char[i,j] <- formatted_val
        } else {
          R_char[i,j] <- ""
        }
      }
    }

    # Paste stars onto correlations
    Rnew <- matrix(paste0(R_char, mystars), ncol = ncol(x_mat))
    diag(Rnew) <- diag(R_char)
    rownames(Rnew) <- colnames(x_mat)

    # Remove triangle
    if (removeTriangle[1] == "upper") {
      Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
    } else {
      Rnew[lower.tri(Rnew, diag = TRUE)] <- ""
    }
    Rnew <- as.data.frame(Rnew, stringsAsFactors = FALSE)

    # Number columns
    colnames(Rnew) <- as.character(seq_len(ncol(Rnew)))
    diag(Rnew)     <- "--"

    # Summary statistics
    if (sumstats) {
      tempstats <- data.frame(
        Mean = colMeans(tempdf, na.rm = TRUE),
        SD   = apply(tempdf, 2, sd, na.rm = TRUE),
        N    = diag(ntemp)
      )
      # Format summary statistics with specified decimals and trailing zeros
      tempstats <- as.data.frame(lapply(tempstats, sprintf, fmt = paste0("%.", decimals, "f")))
      Rnew <- cbind(tempstats, Rnew)
      names(Rnew) <- c("Mean", "SD", "N", seq_len(ncol(R)))
    } else {
      names(Rnew) <- seq_len(ncol(R))
    }

    # Row names: capitalise and number
    rownames(Rnew) <- paste0(
      seq_len(nrow(Rnew)), ". ",
      toupper(substr(rownames(Rnew), 1, 1)),
      substr(rownames(Rnew), 2, nchar(rownames(Rnew)))
    )

    # Attach group label for stacking
    Rnew <- Rnew |>
      tibble::rownames_to_column("Variable") |>
      dplyr::mutate(.group_label = if (!is.null(grp_label)) stringr::str_to_title(grp_label) else "all",
                    .before = Variable)

    Rnew
  }

  # ── Build table data ───────────────────────────────────────────────────────
  if (!is.null(group)) {
    # Split by group, dropping the group column itself
    group_vec  <- x[[group]]
    group_levs <- unique(group_vec)
    x_vars     <- x |> dplyr::select(-dplyr::all_of(group))

    blocks <- purrr::map(group_levs, function(g) {
      df_g <- x_vars[group_vec == g, , drop = FALSE]
      build_block(df_g, grp_label = as.character(g))
    })

    tbl <- dplyr::bind_rows(blocks)
  } else {
    tbl <- build_block(x)
  }

  # ── Return text ────────────────────────────────────────────────────────────
  if (result[1] == "text") {
    return(tbl |> dplyr::select(-.group_label))
  }

  # ── Return LaTeX ───────────────────────────────────────────────────────────
  if (result[1] == "latex") {
    return(xtable::xtable(tbl |> dplyr::select(-.group_label), type = "latex"))
  }

  # ── Build gt table ─────────────────────────────────────────────────────────

  # Row indices per group for tab_row_group
  group_row_map <- tbl |>
    dplyr::mutate(.row = dplyr::row_number()) |>
    dplyr::group_by(.group_label) |>
    dplyr::summarise(start = min(.row), end = max(.row), .groups = "drop")

  out <- tbl |>
    dplyr::select(-.group_label) |>
    gt::gt() |>
    gt::tab_header(title = title) |>
    gt::tab_options(
      table.border.bottom.width = "0px",
      table.border.top.width    = "0px",
      heading.align             = "left"
    ) |>
    gt::cols_align(align = "center", columns = everything()) |>
    gt::cols_align(align = "left",   columns = "Variable") |>
    gt::tab_source_note(
      source_note = gt::html(paste0("<i>Note</i>. ", footer))
    )

  # Add row group labels only when grouping was requested
  if (!is.null(group)) {
    for (i in seq_len(nrow(group_row_map))) {
      out <- out |>
        gt::tab_row_group(
          label = group_row_map$.group_label[i],
          rows  = group_row_map$start[i]:group_row_map$end[i]
        )
    }
  }

  return(out)
}

# Latent Group Model ####
#' This function creates the LGM matrix
#' @param x Data object.
#' @param group Nesting variable.
#' @param title Table capation.
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



