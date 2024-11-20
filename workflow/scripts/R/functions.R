read_config <- function(filename, analysis_name) {
    dat = read_yaml(configfile)[['analyses']]
    found = FALSE
    for (ret in dat) {
        if (ret$analysis == analysis_name) {
            found = TRUE
            break()
        }
    }

    if (found == FALSE) {
        stop('Could not find configuration for analysis "', analysis_name, '"')
    }

    # Consistancy check
    if (is.null(ret$results$names) & is.null(ret$results$contrasts)) {
        stop("At least one value for 'names' or 'contrasts' must be specified under 'results'\n")
    }
    return(ret)
}

read_samplefile = function(filename) {
    ret = fread(filename, check.names=FALSE, data.table=FALSE)

    # Rename samples if requested
    if ("alias" %in% colnames(ret)) {
        rownames(ret) = ret$alias
    } else {
        rownames(ret) = ret$sample
    }

    # Convert strings to factors
    for (col in colnames(ret)) {
        ret[, col] = factor(as.character(ret[, col]), levels=unique(ret[, col]))
    }

    return(ret)
}

read_counts = function(filename, coldat) {
    select = c("character", rep("integer", nrow(coldat)))
    names(select) = c("mirna_id", as.character(coldat$sample))

    ret = fread(
        filename,
        select=select,
        sep=",",
        stringsAsFactors=FALSE,
        check.names=FALSE,
        data.table=FALSE)

    # Add rownames
    rownames(ret) = ret$mirna_id
    ret$mirna_id = NULL

    # Rename samples if requested
    if ("alias" %in% colnames(coldat)) {
        alias = coldat$alias
        names(alias) = coldat$sample
        colnames(ret) = alias[colnames(ret)]
    }

    return(ret)
}

plot_sample_counts = function(dat) {
    cnt = as.data.table(table(dat[, length(which(value >= 1)), by=mirna_id]$V1))
    set(cnt, 1L, 1L, "None")
    cnt[, sample := factor(cnt$V1, levels=cnt$V1)]

    plt = ggplot(cnt, aes(x=sample, y=N, label=N)) +
        geom_bar(stat="identity") +
        geom_text(hjust=-.25, size=2) +
        scale_y_continuous(expand=expansion(mult=c(0, .1))) +
        xlab("Number of samples having a count > 1") +
        ylab("Number of genes") +
        coord_flip() +
        theme_minimal()
    return(plt)
}

plot_group_counts = function(dat) {
    cnt = dat[, length(which(value >= 1)), by=.(mirna_id, group)]
    cnt = melt(cnt[, as.list(table(V1)), by=group], id.var="group")
    cnt[variable == 0, variable := "None"]
    cnt[, variable := factor(cnt$variable, levels=unique(cnt$variable))]

    plt = ggplot(cnt, aes(x=variable, y=value, label=value)) +
        facet_grid(group ~ .) +
        geom_bar(stat="identity") +
        geom_text(hjust=-.25, size=2) +
        scale_y_continuous(expand=expansion(mult=c(0, .1))) +
        xlab("Number of samples having a count > 1") +
        ylab("Number of genes") +
        coord_flip() +
        theme_minimal()

    ggsave("test.pdf", plt, height=4, width=6)
    return(plt)
}

get_strict_counts = function(count, coldat, model) {
    ncounts = 1
    nsample = 4
    ngroups = ifelse("condition" %in% all.vars(model), 4, 2)
    dat = melt(as.data.table(count, keep.rownames="mirna_id"),
        id.vars="mirna_id", variable.name="sample")
    dat = merge(dat, coldat)
    dat[, group := do.call(paste, c(.SD, sep="-")), .SDcols=all.vars(model)]
    p1 = plot_sample_counts(dat)
    p2 = plot_group_counts(dat)
    ggsave("zero-counts.pdf", plot_grid(p1, p2, ncol=1))

    cnt = dat[, length(which(value >= ncounts)), by=.(mirna_id, condition)]
    ids = cnt[V1 >= nsample, .N, by=mirna_id][N >= ngroups, mirna_id] 
    ret = count[ids,]
    return(ret)
}

get_all_vs_all = function(val, coldat) {
    all = combn(levels(coldat[, val[1]]), 2)
    all = apply(all, 2, rev)
    colnames(all) = apply(all, 2, paste0, collapse="-vs-")
    ret = apply(all, 2, function(x) c(val[1], x[1], x[2]), simplify=FALSE)
    return(ret)
}

parse_contrasts = function(dat, dds) {
    coldat = as.data.frame(colData(dds))
    allvars = all.vars(design(dds))

    ret = list()
    for (val in dat) {
        if (val[1] %in% allvars) {
            # Test for c(var, cond, cond)
            if (all(val[2:3] %in% levels(coldat[, val[1]]))) {
                # Both conditions specified
                name = paste0(val[2], "-vs-", val[3])
                ret[[name]] = val
            } else if (val[2] == "all") {
                # Determine all-vs-all comparisons
                ret = c(ret, get_all_vs_all(val, coldat))
            } else {
                # Nothing to do
                stop("Could not determine what comparison to make given:", val, "\n")
            }
        } else if (is.list(val) & length(val) == 1 & !is.null(names(val))) {
            # Have named list
            ret[[names(val)]] = as.list(val[[1]])
        } else {
            stop("Could not determine what comparison to make given:", val, "\n")
        }
    }
    return(ret)
}

parse_names = function(dat, dds) {
    rn = resultsNames(dds)
    ret = list()
    for (val in dat) {
        if (val[2] %in% rn) {
            ret[[val[1]]] = val[2]
        } else {
            stop("Result name \"", val[2], "\" not in resultsNames()\n")
        }
    }
    return(ret)
}

get_comparisons = function(res, dds) {
    contrasts = parse_contrasts(res$contrasts, dds)
    names = parse_names(res$names, dds)
    return(list("contrast"=contrasts, "name"=names))
}

get_results = function(dds, comparisons, config) {
    ret = list()
    cfg = config[setdiff(names(config), c("contrasts", "names"))]
    for (type in names(comparisons)) {
        for (cmp in names(comparisons[[type]])) {
            val = list(comparisons[[type]][[cmp]])
            names(val) = type
            res = do.call(results, c(list(dds), val, cfg))
            ret[[cmp]] = lfcShrink(dds, res=res, type="ashr")
        }
    }
    return(ret)
}

get_group_counts = function(dds, config) {
    if (!is.null(config$group_counts)) {
        groups = config$group_counts
    } else {
        groups = all.vars(design(dds))
    }
    cnt = counts(dds, normalized=TRUE)
    dat = as.data.frame(colData(dds))
    ret_list = list()
    for (grp in groups) {
        for (lvl in levels(dat[, grp])) {
            ids = which(dat[, grp] == lvl)
            if (length(ids) > 1) {
                ret_list[[lvl]] = rowMeans(cnt[, ids])
            } else {
                ret_list[[lvl]] = count[, ids]
            }
        }
    }
    ret = do.call(data.frame, ret_list)
    names(ret) = paste0(names(ret_list), "_mean")
    ret$mirna_id = rownames(ret)
    return(ret)
}

format_tab = function(res, dds, cnt, genes=NA) {
    dat = as.data.frame(colData(dds))
    res = as.data.frame(res)
    res$mirna_id = rownames(res)
    ret = data.table("mirna_id"=rownames(res))

    #ret = merge(ret, genes, all.x=TRUE)
    ret = merge(ret, res)
    ret = merge(ret, cnt)
    ret[, baseMean := NULL]
    ret = ret[order(ret$pvalue, na.last=TRUE),]
    return(ret)
}

write_xls = function(tab, cmps, filename) {
    ret = list()
    names = cmps$name
    contrasts = cmps$contrast

    # Look for contrasts / coefs
    for (cmp in names(tab)) {
      ret[[cmp]] = copy(tab[[cmp]])

      if (cmp %in% names) {
        lfc_name = paste0("logFC (", cmp, ")")
      } else {
        vals = contrasts[[cmp]]
        if (!is.null(vals)) {
          if (is.list(vals)) {
            lfc_name = paste0("logFC (", cmp, ")")
          } else {
            lfc_name = paste0("logFC (", vals[2], " / ", vals[3], ")")
          }
        }
      }
      setnames(ret[[cmp]], "log2FoldChange", lfc_name)
    }

    desc_vec = c(
        "mirna_id"="miRNA identifier",
        "logFC (group-A / group-B)"="Shrunken estimates of log2 fold change between two groups (effect size, see: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#moreshrink)",
        "lfcSE"="Standard error of logFC column",
        "pvalue"="p-value of test for features passing independent filtering (see: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#indfilt)",
        "padj"="p-value adjusted for multiple testing using method proposed by Benjamini, Y., and Hochberg, Y. (1995).",
        "<group-A>_mean"="Mean of normalized counts among samples associated with group-A",
        "<group-B>_mean"="Mean of normalized counts among samples associated with group-B"
    )

    desc = data.table("Columns"=names(desc_vec), "Description"=desc_vec)
    ret = c(list("Column Descriptions"=desc), ret)
    WriteXLS(ret, ExcelFileName=filename, AdjWidth=TRUE, BoldHeaderRow=TRUE, FreezeRow=1)
}

write_counts = function(mat, filename) {
    ret = data.table("mirna_id"=rownames(mat), mat)
    fwrite(ret, file=filename)
}

# Plotting
plot_pca = function(
    mat, dat, pcs=c("PC1", "PC2"),
    color_by=NULL, shape_by=NULL, 
    label_by=NULL, size=NULL, alpha=NULL,
    filename=NULL) {

    require(ggplot2)
    require(ggrepel)
    require(data.table)

    # Set defaults
    size <- ifelse(is.null(size), 3.5, size)
    alpha <- ifelse(is.null(alpha), 1, alpha)

    # rlog values
    mat = assay(rlog(dds))

    # Compute pca and add to colData
    pca <- prcomp(t(mat))
    pca_x <- data.table(pca$x, keep.rownames=TRUE)
    dat <- data.table(dat, keep.rownames=TRUE)
    plt <- merge(pca_x, dat)

    # Get percent variance
    pctvar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 2)
    pctvar <- pctvar[which(colnames(pca$x) %in% pcs)]

    # Make base plot
    ret <- ggplot(plt, aes_string(x=pcs[1], y=pcs[2], color=color_by, shape=shape_by, label=label_by))
    ret <- ret + xlab(paste0(pcs[1], " (", pctvar[1], "% variance)"))
    ret <- ret + ylab(paste0(pcs[2], " (", pctvar[2], "% variance)"))
    ret <- ret + geom_point(size=size, alpha=alpha)
    ret <- ret + theme_light()

    if (!is.null(label_by)) {
        ret <- ret + geom_text_repel(
            aes_string(label=label_by),
            size=2.75,
            box.padding=0.5,
            point.padding=0.5,
            min.segment.length=0,
            max.overlaps=Inf,
            show.legend=FALSE)
    }

    # Save if requested
    if (!is.null(filename)) {
        ggsave(filename=filename, plot=ret, width=7.5, height=6)
    }

    # Return plot if requested
    invisible(ret)
}

get_contrast_samples = function(dat, cmp, comparisons) {
    ctr = comparisons$contrast[[cmp]]
    if (!is.null(ctr) & ctr[1] %in% colnames(dat)) {
        ret = rownames(dat[which(dat[,ctr[1]] %in% ctr[2:3]),, drop=TRUE])
    } else {
        ret = rownames(dat)
    }
    return(ret)
}

is_true = function(x) {
    ret = ifelse(!is.null(x) & x == TRUE, TRUE, FALSE)
    return(ret)
}

plot_heatmap = function(res, cmp, dds, config, comparisons, filename, n=100) {
    # Data to plot
    dat = as.data.frame(colData(dds))[, all.vars(design(dds)), drop=FALSE]
    res = na.omit(as.data.frame(res[[cmp]]))
    res = res[order(abs(res$log2FoldChange), decreasing=TRUE, na.last=TRUE),]

    # Subset ids to use (top 100, ranked by LFC)
    ids = na.omit(rownames(res)[1:100])

    # Determine which samples to plot
    if (is_true(config$heatmap$only_contrasts)) {
        samples = get_contrast_samples(dat, cmp, comparisons)
        dat = dat[samples,]
    } else {
        samples = rownames(dat)
    }

    # Z-score rlog counts
    mat = assay(rlog(dds))[ids, samples]
    mat = t(scale(t(mat)))

    # zlim
    lim = max(abs(range(mat)))
    breaks = seq(-lim, lim, length.out=101)

    # Use spectral but replace midpoint with white
    colors = brewer.pal(11, "Spectral")
    colors[6] = "#FFFFFF"
    colors = colorRampPalette(colors)(101)

    # Plot
    pheatmap(
        mat, colors, annotation_col=dat, breaks=breaks, scale="none", border_color=NA,
        show_rownames=FALSE, show_colnames=TRUE, angle_col=90,
        cluster_rows=TRUE, clustering_distance_rows="correlation",
        cluster_cols=TRUE, clustering_distance_cols="correlation",
        treeheight_row=0, treeheight_col=30,
        cellwidth=20, filename=filename, legend_col="")
}

plot_dispersions = function(dds, filename) {
    pdf(file=filename, height=6, width=8)
    plotDispEsts(dds)
    dev.off()
}

plot_sparsity = function(dds, filename) {
    pdf(file=filename, height=6, width=8)
    plotSparsity(dds)
    dev.off()
}

plot_volcano = function(res, filename, ntop=25, genes=NA) {
    require(ggrepel)
    # Ablines
    vline = c(-2, 2)
    hline = -log10(0.05)
    res = as.data.frame(copy(res))
    dat = as.data.table(res, keep.rownames="mirna_id")

    # Filter NA
    dat = dat[
        !is.na(pvalue) &
        !is.na(padj) &
        !is.na(log2FoldChange)]

    # Format symbol
    dat[, symbol := sub("^.*-", "", mirna_id)]

    # Add stats and do selection
    dat[, "pval" := -log10(pvalue)]
    dat[, "alfc" := abs(log2FoldChange)]
    setorder(dat, padj, -alfc)
    dat[, "rank" := nrow(dat):1]
    dat[!1:ifelse(nrow(dat) < ntop, nrow(dat), ntop), "symbol" := NA]

    # Plot
    ret = ggplot(dat, aes(x=log2FoldChange, y=pval)) +
        geom_vline(xintercept=vline, linetype="dashed", size=0.3, color="grey20") +
        geom_hline(yintercept=hline, linetype="dashed", size=0.3, color="grey20") +
        geom_point(alpha=0.3, size=1.5, shape=16) +
        geom_point(data=dat[1:ntop,], size=1.5, shape=16, color="orangered") +
        geom_label_repel(
            data=dat[1:ntop,],
            aes(label=symbol),
            size=2.5,
            max.overlaps=Inf,
            min.segment.length=0,
            point.padding=0.25,
            box.padding=0.6,
            label.padding=0.25,
            force_pull=0,
            max.time=5,
            max.iter=1e6,
            na.rm=TRUE,
            segment.size=0.05,
            verbose=TRUE,
            show.legend=FALSE) +
        xlab("log2 Fold-Change") +
        ylab("-log10(p-value)") +
        theme_minimal()
    ggsave(ret, file=filename, width=7, height=6)
}

