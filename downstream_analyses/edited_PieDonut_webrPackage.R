function (data, mapping, start = getOption("PieDonut.start", 
    0), addPieLabel = TRUE, addDonutLabel = TRUE, showRatioDonut = TRUE, 
    showRatioPie = TRUE, ratioByGroup = TRUE, showRatioThreshold = getOption("PieDonut.showRatioThreshold", 
        0.02), labelposition = getOption("PieDonut.labelposition", 
        2), labelpositionThreshold = 0.1, r0 = getOption("PieDonut.r0", 
        0.3), r1 = getOption("PieDonut.r1", 1), r2 = getOption("PieDonut.r2", 
        1.2), explode = NULL, selected = NULL, explodePos = 0.1, 
    color = "white", pieAlpha = 0.8, donutAlpha = 1, maxx = NULL, 
    showPieName = TRUE, showDonutName = FALSE, title = NULL, 
    pieLabelSize = 4, donutLabelSize = 3, titlesize = 5, explodePie = TRUE, 
    explodeDonut = FALSE, use.label = TRUE, use.labels = TRUE, 
    family = getOption("PieDonut.family", "")) 
{
    (cols = colnames(data))
    if (use.labels) 
        data = addLabelDf(data, mapping)
    count <- NULL
    if ("count" %in% names(mapping)) 
        count <- getMapping(mapping, "count")
    count
    pies <- donuts <- NULL
    (pies = getMapping(mapping, "pies"))
    if (is.null(pies)) 
        (pies = getMapping(mapping, "pie"))
    if (is.null(pies)) 
        (pies = getMapping(mapping, "x"))
    (donuts = getMapping(mapping, "donuts"))
    if (is.null(donuts)) 
        (donuts = getMapping(mapping, "donut"))
    if (is.null(donuts)) 
        (donuts = getMapping(mapping, "y"))
    if (!is.null(count)) {
        df <- data %>% group_by(.data[[pies]]) %>% dplyr::summarize(Freq = sum(.data[[count]]))
        df
    }
    else {
        df = data.frame(table(data[[pies]]))
    }
    colnames(df)[1] = pies
    df$end = cumsum(df$Freq)
    df$start = dplyr::lag(df$end)
    df$start[1] = 0
    total = sum(df$Freq)
    df$start1 = df$start * 2 * pi/total
    df$end1 = df$end * 2 * pi/total
    df$start1 = df$start1 + start
    df$end1 = df$end1 + start
    df$focus = 0
    if (explodePie) 
        df$focus[explode] = explodePos
    df$mid = (df$start1 + df$end1)/2
    df$x = ifelse(df$focus == 0, 0, df$focus * sin(df$mid))
    df$y = ifelse(df$focus == 0, 0, df$focus * cos(df$mid))
    df$label = df[[pies]]
    df$ratio = df$Freq/sum(df$Freq)
    if (showRatioPie) {
        df$label = ifelse(df$ratio >= showRatioThreshold, paste0(df$label, 
            "\n(", scales::percent(df$ratio), ")"), as.character(df$label))
    }
    df$labelx = (r0 + r1)/2 * sin(df$mid) + df$x
    df$labely = (r0 + r1)/2 * cos(df$mid) + df$y
    if (!is.factor(df[[pies]])) 
        df[[pies]] <- factor(df[[pies]])
    df
    mainCol = c("#e6550d", "#31a354", "#de2d26", "#3182bd", "#2c7fb8")
    df$radius = r1
    df$radius[df$focus != 0] = df$radius[df$focus != 0] + df$focus[df$focus != 
        0]
    df$hjust = ifelse((df$mid%%(2 * pi)) > pi, 1, 0)
    df$vjust = ifelse(((df$mid%%(2 * pi)) < (pi/2)) | (df$mid%%(2 * 
        pi) > (pi * 3/2)), 0, 1)
    df$segx = df$radius * sin(df$mid)
    df$segy = df$radius * cos(df$mid)
    df$segxend = (df$radius + 0.05) * sin(df$mid)
    df$segyend = (df$radius + 0.05) * cos(df$mid)
    df
    if (!is.null(donuts)) {
        subColor = makeSubColor(mainCol, no = length(unique(data[[donuts]])))
        subColor
        data
        if (!is.null(count)) {
            df3 <- as.data.frame(data[c(donuts, pies, count)])
            colnames(df3) = c("donut", "pie", "Freq")
            df3
            df3 <- eval(parse(text = "complete(df3,donut,pie)"))
            df3$Freq[is.na(df3$Freq)] = 0
            if (!is.factor(df3[[1]])) 
                df3[[1]] = factor(df3[[1]])
            if (!is.factor(df3[[2]])) 
                df3[[2]] = factor(df3[[2]])
            df3 <- df3 %>% arrange(.data$pie, .data$donut)
            a <- df3 %>% spread(.data$pie, value = .data$Freq)
            a = as.data.frame(a)
            a
            rownames(a) = a[[1]]
            a = a[-1]
            a
            colnames(df3)[1:2] = c(donuts, pies)
        }
        else {
            df3 = data.frame(table(data[[donuts]], data[[pies]]), 
                stringsAsFactors = FALSE)
            colnames(df3)[1:2] = c(donuts, pies)
            a = table(data[[donuts]], data[[pies]])
            a
        }
        a
        df3
        df3$group = rep(colSums(a), each = nrow(a))
        df3$pie = rep(1:ncol(a), each = nrow(a))
        total = sum(df3$Freq)
        total
        df3$ratio1 = df3$Freq/total
        df3
        if (ratioByGroup) {
            df3$ratio = scales::percent(df3$Freq/df3$group)
        }
        else {
            df3$ratio <- scales::percent(df3$ratio1)
        }
        df3$end = cumsum(df3$Freq)
        df3
        df3$start = dplyr::lag(df3$end)
        df3$start[1] = 0
        df3$start1 = df3$start * 2 * pi/total
        df3$end1 = df3$end * 2 * pi/total
        df3$start1 = df3$start1 + start
        df3$end1 = df3$end1 + start
        df3$mid = (df3$start1 + df3$end1)/2
        df3$focus = 0
        if (!is.null(selected)) {
            df3$focus[selected] = explodePos
        }
        else if (!is.null(explode)) {
            selected = c()
            for (i in 1:length(explode)) {
                start = 1 + nrow(a) * (explode[i] - 1)
                selected = c(selected, start:(start + nrow(a) - 
                  1))
            }
            selected
            df3$focus[selected] = explodePos
        }
        df3
        df3$x = 0
        df3$y = 0
        df
        if (!is.null(explode)) {
            explode
            for (i in 1:length(explode)) {
                xpos = df$focus[explode[i]] * sin(df$mid[explode[i]])
                ypos = df$focus[explode[i]] * cos(df$mid[explode[i]])
                df3$x[df3$pie == explode[i]] = xpos
                df3$y[df3$pie == explode[i]] = ypos
            }
        }
        df3$no = 1:nrow(df3)
        df3$label = df3[[donuts]]
        if (showRatioDonut) {
            if (max(nchar(levels(df3$label))) <= 2) 
                df3$label = paste0(df3$label, "(", df3$ratio, 
                  ")")
            else df3$label = paste0(df3$label, "\n(", df3$ratio, 
                ")")
        }
        df3$label[df3$ratio1 == 0] = ""
        df3$label[df3$ratio1 < showRatioThreshold] = ""
        df3$hjust = ifelse((df3$mid%%(2 * pi)) > pi, 1, 0)
        df3$vjust = ifelse(((df3$mid%%(2 * pi)) < (pi/2)) | (df3$mid%%(2 * 
            pi) > (pi * 3/2)), 0, 1)
        df3$no = factor(df3$no)
        df3
        labelposition
        if (labelposition > 0) {
            df3$radius = r2
            if (explodeDonut) 
                df3$radius[df3$focus != 0] = df3$radius[df3$focus != 
                  0] + df3$focus[df3$focus != 0]
            df3$segx = df3$radius * sin(df3$mid) + df3$x
            df3$segy = df3$radius * cos(df3$mid) + df3$y
            df3$segxend = (df3$radius + 0.05) * sin(df3$mid) + 
                df3$x
            df3$segyend = (df3$radius + 0.05) * cos(df3$mid) + 
                df3$y
            if (labelposition == 2) 
                df3$radius = (r1 + r2)/2
            df3$labelx = (df3$radius) * sin(df3$mid) + df3$x
            df3$labely = (df3$radius) * cos(df3$mid) + df3$y
        }
        else {
            df3$radius = (r1 + r2)/2
            if (explodeDonut) 
                df3$radius[df3$focus != 0] = df3$radius[df3$focus != 
                  0] + df3$focus[df3$focus != 0]
            df3$labelx = df3$radius * sin(df3$mid) + df3$x
            df3$labely = df3$radius * cos(df3$mid) + df3$y
        }
        df3$segx[df3$ratio1 == 0] = 0
        df3$segxend[df3$ratio1 == 0] = 0
        df3$segy[df3$ratio1 == 0] = 0
        df3$segyend[df3$ratio1 == 0] = 0
        if (labelposition == 0) {
            df3$segx[df3$ratio1 < showRatioThreshold] = 0
            df3$segxend[df3$ratio1 < showRatioThreshold] = 0
            df3$segy[df3$ratio1 < showRatioThreshold] = 0
            df3$segyend[df3$ratio1 < showRatioThreshold] = 0
        }
        df3
        del = which(df3$Freq == 0)
        del
        if (length(del) > 0) 
            subColor <- subColor[-del]
        subColor
    }
    p <- ggplot() + theme_no_axes() + coord_fixed()
    if (is.null(maxx)) {
        r3 = r2 + 0.3
    }
    else {
        r3 = maxx
    }
    p1 <- p + geom_arc_bar(aes_string(x0 = "x", y0 = "y", r0 = as.character(r0), 
        r = as.character(r1), start = "start1", end = "end1", 
        fill = pies), alpha = pieAlpha, color = color, data = df) + 
        transparent() + scale_fill_manual(values = mainCol) + 
        xlim(r3 * c(-1, 1)) + ylim(r3 * c(-1, 1)) + guides(fill = FALSE)
    if ((labelposition == 1) & (is.null(donuts))) {
        p1 <- p1 + geom_segment(aes_string(x = "segx", y = "segy", 
            xend = "segxend", yend = "segyend"), data = df) + 
            geom_text(aes_string(x = "segxend", y = "segyend", 
                label = "label", hjust = "hjust", vjust = "vjust"), 
                color = "black", size = pieLabelSize, data = df, 
                family = family)
    }
    else if ((labelposition == 2) & (is.null(donuts))) {
        p1 <- p1 + geom_segment(aes_string(x = "segx", y = "segy", 
            xend = "segxend", yend = "segyend"), data = df[df$ratio < 
            labelpositionThreshold, ]) + geom_text(aes_string(x = "segxend", 
            y = "segyend", label = "label", hjust = "hjust", 
            vjust = "vjust"), color = "black", size = pieLabelSize, 
            data = df[df$ratio < labelpositionThreshold, ], family = family) + 
            geom_text(aes_string(x = "labelx", y = "labely", 
                label = "label"), color = "black", size = pieLabelSize, 
                data = df[df$ratio >= labelpositionThreshold, 
                  ], family = family)
    }
    else {
        p1 <- p1 + geom_text(aes_string(x = "labelx", y = "labely", 
            label = "label"), color = "black", size = pieLabelSize, 
            data = df, family = family)
    }
    if (showPieName) 
        p1 <- p1 + annotate("text", x = 0, y = 0, label = pies, 
            size = titlesize, family = family)
    p1 <- p1 + theme(text = element_text(family = family))
    if (!is.null(donuts)) {
        if (explodeDonut) {
            p3 <- p + geom_arc_bar(aes_string(x0 = "x", y0 = "y", 
                r0 = as.character(r1), r = as.character(r2), 
                start = "start1", end = "end1", fill = "no", 
                explode = "focus"), alpha = donutAlpha, color = color, 
                data = df3)
        }
        else {
            p3 <- p + geom_arc_bar(aes_string(x0 = "x", y0 = "y", 
                r0 = as.character(r1), r = as.character(r2), 
                start = "start1", end = "end1", fill = "no"), 
                alpha = donutAlpha, color = color, data = df3)
        }
        p3 <- p3 + transparent() + scale_fill_manual(values = subColor) + 
            xlim(r3 * c(-1, 1)) + ylim(r3 * c(-1, 1)) + guides(fill = FALSE)
        p3
        if (labelposition == 1) {
            p3 <- p3 + geom_segment(aes_string(x = "segx", y = "segy", 
                xend = "segxend", yend = "segyend"), data = df3) + 
                geom_text(aes_string(x = "segxend", y = "segyend", 
                  label = "label", hjust = "hjust", vjust = "vjust"), 
                  color = "black", size = donutLabelSize, data = df3, 
                  family = family)
        }
        else if (labelposition == 0) {
            p3 <- p3 + geom_text(aes_string(x = "labelx", y = "labely", 
                label = "label"), color = "black", size = donutLabelSize, 
                data = df3, family = family)
        }
        else {
            p3 <- p3 + geom_segment(aes_string(x = "segx", y = "segy", 
                xend = "segxend", yend = "segyend"), data = df3[df3$ratio1 < 
                labelpositionThreshold, ]) + geom_text(aes_string(x = "segxend", 
                y = "segyend", label = "label", hjust = "hjust", 
                vjust = "vjust"), color = "black", size = donutLabelSize, 
                data = df3[df3$ratio1 < labelpositionThreshold, 
                  ], family = family) + geom_text(aes_string(x = "labelx", 
                y = "labely", label = "label"), color = "black", 
                size = donutLabelSize, data = df3[df3$ratio1 >= 
                  labelpositionThreshold, ], family = family)
        }
        if (!is.null(title)) 
            p3 <- p3 + annotate("text", x = 0, y = r3, label = title, 
                size = titlesize, family = family, color = "black")
        else if (showDonutName) 
            p3 <- p3 + annotate("text", x = (-1) * r3, y = r3, 
                label = donuts, hjust = 0, size = titlesize, 
                family = family, color = "black")
        p3 <- p3 + theme(text = element_text(family = family, 
            color = "black"))
        grid.newpage()
        p1 <- p1 + theme(panel.background = element_rect(fill = "transparent", 
            colour = NA), plot.background = element_rect(fill = "transparent", 
            colour = NA), panel.grid.minor = element_blank(), 
            panel.grid.major = element_blank())
        p3 <- p3 + theme(panel.background = element_rect(fill = "transparent", 
            colour = NA), plot.background = element_rect(fill = "transparent", 
            colour = NA), panel.grid.minor = element_blank(), 
            panel.grid.major = element_blank())
        print(p3, vp = viewport(height = 1, width = 1))
        print(p1, vp = viewport(height = 1, width = 1))
    }
    else {
        p1
    }
}
