# This script is written for analyzing a dyadic conversation (i.e., between therapist and client).
# In addition, note that it assumes a therpaist turn comes first.

# setup -------------------------------------------------------------------

library(tidyverse)
library(stringr)
library(tidytext) # 'tm' package also needs to be installed for 'cast_tdm' function to run
library(textstem)
library(reshape2)

# some relevant information retrieved from 'sessionInfo()'
# R version 3.6.0 (2019-04-26)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 17134)
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] reshape2_1.4.3       textstem_0.1.4       koRpus.lang.en_0.1-2 koRpus_0.11-5        sylly_0.1-5         
# [6] tidytext_0.2.0       forcats_0.4.0        stringr_1.4.0        dplyr_0.8.0.1        purrr_0.3.2         
# [11] readr_1.3.1          tidyr_0.8.3          tibble_2.1.1         ggplot2_3.1.1        tidyverse_1.2.1  

# function ----------------------------------------------------------------

get_turn_sim <- function(v) {
  # v = feature matrix
  s <- matrix(0, nrow = ncol(v), ncol = ncol(v))
  for (i in 1:ncol(v)) {
    for (j in 1:i) {
        s[i, j] <- (v[, i] %*% v[, j]) / (sqrt(sum(v[, i]^2)*sum(v[, j]^2)))
      }
      # message("word index in progress: ", "(", i, ", ", j, ")")
  }
  s[is.nan(s)] <- 0
  s
}

get_word_sim <- function(o, c, w) {
  # o = occurence vector, c = cooccurence matrix, w = word list with sentence window
  n <- max(w[2])
  s <- matrix(nrow =  nrow(c), ncol = ncol(c), dimnames = dimnames(c))
  oo <- map(seq_along(o[[1]]), function(x) o[[2]][x] + o[[2]]) # make a list for all possible combinations of sum
  oo <- matrix(unlist(oo, use.names = FALSE), ncol = nrow(o)) # transform a list into a matrix
  cond_1 <- oo == c + n
  cond_2 <- o[[2]] == c
  p_no_bar <- c / n
  p_both_bar <- ifelse (cond_1, 1, (n - oo + c) / n)
  p_left_bar <- ifelse (cond_2, 1, (o[[2]] - c) / n)
  p_right_bar <- ifelse (t(cond_2), 1, t((o[[2]] - c) / n))
  s <- (p_no_bar * p_both_bar) / (p_left_bar * p_right_bar)
  # According to the similarity algorithm presented by Angus, Smith, and Wiles (2012),
  # similarity values between the same term are not calculated properly
  # (i.e., they are smaller than other similarity values between different terms).
  # Therefore, in this application, similarity values between the same terms are replaced by
  # the maximum similarity value between different terms.
  diag(s) <- max(s)
  s
}

get_feature_matrix <- function(data) {
  # list unique words in each sentence window
  unique_words_win <- data %>%
    select(word, win_num) %>% 
    distinct()
  
  # count the frequency
  unique_wordcount_win <- unique_words_win %>% 
    count(word) # 'sort = TRUE' if you want to see the top words
  
  # create co-occurence matrix  (https://stackoverflow.com/questions/13281303/creating-co-occurrence-matrix)
  co_mat <- crossprod(table(unique_words_win[2:1]))
  
  word_sim <- get_word_sim(unique_wordcount_win, co_mat, unique_words_win)
  
  # extract keywords // threshold to be keywords can be flexible (In this script, it's 50).
  keyword_50 <- data %>% 
    count(word) %>% 
    top_n(50) %>% 
    .$word
  
  # similarity matrix
  sim_matrix <- word_sim[keyword_50, ]
  
  # boolean matrix
  turn_word <- data %>%
    select(word, turn) %>% 
    distinct()
  boolean_matrix_wo <- turn_word %>% 
    group_by(turn) %>% 
    count(word) %>% 
    cast_tdm(word, turn, n)
  boolean_matrix_wo <- as.matrix(boolean_matrix_wo)
  boolean_matrix_wo <- boolean_matrix_wo[order(rownames(boolean_matrix_wo)), ]
  # add columns with no words
  # from https://stackoverflow.com/questions/42978822/r-add-missing-columns-and-rows-of-data-dplyr-tidyr-complete
  boolean_matrix <- matrix(0, nrow = nrow(unique_wordcount_win), ncol = nrow(data_raw), dimnames = list(rownames(boolean_matrix_wo), 1:nrow(data_raw)))
  boolean_matrix[row.names(boolean_matrix_wo), colnames(boolean_matrix_wo)] <- boolean_matrix_wo
  
  # feature matrix
  feature_matrix <- sim_matrix %*% boolean_matrix
}

get_crp <- function(data) {
  # generate 'data_tc' to visualize triangles for 'Therapist-Client' interactions.
  data_tc <- data %>% 
    filter(category == "TC") %>% 
    mutate(x1 = x.posit - (w.size/2), 
           x2 = x.posit + (w.size/2),
           y1 = y.posit - (h.size/2),
           y2 = y.posit + (h.size/2)) %>%
    unite(x1, y1, col = "a", remove = FALSE) %>% 
    unite(x1, y2, col = "b", remove = FALSE) %>% 
    unite(x2, y1, col = "c", remove = FALSE) %>% 
    unite(x2, y2, col = "d", remove = FALSE) %>%
    mutate(e = b, f = c) %>% 
    gather(a, b, c, d, e, f, key = "key", value = "temp") %>% 
    separate(temp, into = c("x", "y")) %>%
    mutate(x = as.numeric(x), y = as.numeric(y)) %>% 
    arrange(Var2, Var1)
  data_tc$group_id <- data_tc %>% group_indices(Var2, Var1)
  data_tc <- data_tc %>% 
    group_by(group_id) %>% 
    mutate(split_half = if_else(row_number() <= 3, 1, 2)) %>% 
    mutate(category = if_else(split_half == 1, "T", 
                              if_else(split_half == 2, "C", as.character(category)))) %>%
    ungroup() %>% 
    mutate(category = as_factor(category)) %>% 
    select(-(x.posit:key))
  data_tc$group_id <- data_tc %>% 
    group_indices(group_id, split_half)
  # generate CRP
  # combine 'geom_tile' (for squares (T-T & C-C interactions)) and 'geom_polygon' (for triangles (T-C interactions))
  data %>% 
      filter(category != "TC") %>% 
      ggplot() + 
      geom_tile(aes(x = x.posit, y = y.posit, alpha = value, fill = category, width = w.size, height = h.size)) +
      geom_polygon(data = data_tc, aes(x = x, y = y, alpha = value, fill = category, group = group_id)) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_reverse(expand = c(0, 0)) + 
      scale_alpha(range = c(0, 1)) +
      scale_fill_manual(values = c("T" = "#0072B2", "C" = "#D55E00", "N" = "#ffffff"), 
                        breaks = c("T", "C", "N"), labels = c("Therapist", "Client", "")) +
      labs(x = "", y = "") +
      guides(alpha = FALSE) + 
      theme_classic(base_size = 11) + 
      theme(aspect.ratio = 1, legend.title = element_blank(), 
            axis.text = element_blank(), axis.ticks = element_blank(),legend.position = c(0.9, 0.9))
}

# calculate primatives for MPR metrics
# written only for a dyadic conversation
# will only generate a normalized version of primitives
get_prim <- function(data, scale, direction, type, participant, medium_value = 5) {
  if (!is.matrix(data)) stop ("Please check the data class. It must be a matrix.")
  filter_list <- list(seq(1, ncol(data), 2), seq(2, ncol(data), 2))
  part_filter <- if ((participant == 1) | (participant == 2)) {
    filter_list[[participant]]
  } else {
    stop ("You put an invalid participant name")
  }
  type_filter <- if (type == "self") {
    filter_list[[participant]]
  } else if (type == "other") {
    filter_list[[-participant]]
  } else {
    stop ("You put an invalid type name")
  }
  s_length <- if (scale == "short") {
    1
  } else if (scale == "medium") {
    medium_value
  } else if (scale == "long") {
    max(map_int(filter_list, length))
  } else {
    stop ("You put an invalid scale name")
  }
  output <- numeric(length = ncol(data))
  for(i in part_filter) {
    if (direction == "forward") {
      output[i] <- sum(data[head(type_filter[type_filter > i], s_length), i]) / 
        length(data[head(type_filter[type_filter > i], s_length), i])
    } else if (direction == "backward") {
      output[i] <- sum(data[i, tail(type_filter[type_filter < i], s_length)]) / 
        length(data[i, tail(type_filter[type_filter < i], s_length)])
    } else {
      stop ("You put an invalid direction name")
    }
  }
  output <- if (participant == 1) {
    output[c(TRUE, FALSE)]
  } else {
    output[!c(TRUE, FALSE)]
  }
  output
}

# return a list of mpr metrics
get_mpr_metrics <- function(data, participant) {
  sbs <- get_prim(data, "short", "backward", "self", participant)
  sfs <- get_prim(data, "short", "forward", "self", participant)
  sbo <- get_prim(data, "short", "backward", "other", participant)
  sfo <- get_prim(data, "short", "forward", "other", participant)
  mbs <- get_prim(data, "medium", "backward", "self", participant)
  mfs <- get_prim(data, "medium", "forward", "self", participant)
  mbo <- get_prim(data, "medium", "backward", "other", participant)
  mfo <- get_prim(data, "medium", "forward", "other", participant)
  lbs <- get_prim(data, "long", "backward", "self", participant)
  lfs <- get_prim(data, "long", "forward", "self", participant)
  lbo <- get_prim(data, "long", "backward", "other", participant)
  lfo <- get_prim(data, "long", "forward", "other", participant)
  metrics <- tibble(sbs = sbs, sfs = sfs, sbo = sbo, sfo = sfo, 
                        mbs = mbs, mfs = mfs, mbo = mbo, mfo = mfo, 
                        lbs = lbs, lfs = lfs, lbo = lbo, lfo = lfo)
  metrics[is.na(metrics)] <- 0
  metrics <- metrics %>% mutate(participant = participant, 
                                ITR = sbo, TI = mfo * (1 - sbo), TR = mbo * sbo, TCO = mbo + mfo,
                                TCS = mbs + mfs, LTN = lfo - lbo, LTCO = lfo + lbo, LTCS = lfs + lbs)
}

# plot changes of mpr metrics during the conversation
get_mpr_plot <- function(data, mpr) {
  ggplot(data, aes_string(x = "turn", y = deparse(substitute(mpr)))) +
    geom_line(aes(color = participant, linetype=participant)) + 
    geom_smooth(aes(color = participant, linetype=participant), se = FALSE, method = "loess") +
    scale_color_manual(values = c("1" = "#0072B2", "2" = "#D55E00"), 
                       labels = c("Therapist", "Client")) +
    scale_linetype_manual(values=c("1" = "solid", "2" = "longdash"), 
                          labels = c("Therapist", "Client"))
}

# contractions ------------------------------------------------------------
#Adapted from 'https://stackoverflow.com/questions/19790188/expanding-english-language-contractions-in-python'
contractions_list <- 
  c("ain't" = "am not",
    "aren't" = "are not",
    "can't" = "cannot",
    "can't've" = "cannot have",
    "'cause" = "because",
    "could've" = "could have",
    "couldn't" = "could not",
    "couldn't've" = "could not have",
    "didn't" = "did not",
    "doesn't" = "does not",
    "don't" = "do not",
    "hadn't" = "had not",
    "hadn't've" = "had not have",
    "hasn't" = "has not",
    "haven't" = "have not",
    "he'd" = "he would", #"he had / he would"
    "he'd've" = "he would have",
    "he'll" = "he will", #"he shall / he will"
    "he'll've" = "he will have", #"he shall have / he will have"
    "he's" = "he is", #"he has / he is"
    "how'd" = "how did",
    "how'd'y" = "how do you",
    "how'll" = "how will",
    "how's" = "how is", #"how has / how is / how does"
    "i'd" = "i would", #"I had / I would"
    "i'd've" = "i would have",
    "i'll" = "i will", #"I shall / I will"
    "i'll've" = "i will have", #"I shall have / I will have"
    "i'm" = "i am",
    "i've" = "i have",
    "isn't" = "is not",
    "it'd" = "it would", #"it had / it would"
    "it'd've" = "it would have",
    "it'll" = "it will", #"it shall / it will"
    "it'll've" = "it will have", #"it shall have / it will have"
    "it's" = "it is", #"it has / it is"
    "let's" = "let us",
    "ma'am" = "madam",
    "mayn't" = "may not",
    "might've" = "might have",
    "mightn't" = "might not",
    "mightn't've" = "might not have",
    "must've" = "must have",
    "mustn't" = "must not",
    "mustn't've" = "must not have",
    "needn't" = "need not",
    "needn't've" = "need not have",
    "o'clock" = "of the clock",
    "oughtn't" = "ought not",
    "oughtn't've" = "ought not have",
    "shan't" = "shall not",
    "sha'n't" = "shall not",
    "shan't've" = "shall not have",
    "she'd" = "she would", #"she had / she would"
    "she'd've" = "she would have",
    "she'll" = "she will", #"she shall / she will"
    "she'll've" = "she will have", #"she shall have / she will have"
    "she's" = "she is", #"she has / she is"
    "should've" = "should have",
    "shouldn't" = "should not",
    "shouldn't've" = "should not have",
    "so've" = "so have",
    "so's" = "so is", #"so as / so is"
    "that'd" = "that would", #"that would / that had"
    "that'd've" = "that would have",
    "that's" = "that is", #"that has / that is"
    "there'd" = "there would", #"there had / there would"
    "there'd've" = "there would have",
    "there's" = "there is", #"there has / there is"
    "they'd" = "they would", #"they had / they would"
    "they'd've" = "they would have",
    "they'll" = "they will", #"they shall / they will"
    "they'll've" = "they will have", #"they shall have / they will have"
    "they're" = "they are",
    "they've" = "they have",
    "to've" = "to have",
    "wasn't" = "was not",
    "we'd" = "we would", #"we had / we would"
    "we'd've" = "we would have",
    "we'll" = "we will",
    "we'll've" = "we will have",
    "we're" = "we are",
    "we've" = "we have",
    "weren't" = "were not",
    "what'll" = "what will", #"what shall / what will"
    "what'll've" = "what will have", #"what shall have / what will have"
    "what're" = "what are",
    "what's" = "what is", #"what has / what is"
    "what've" = "what have",
    "when's" = "when is", #"when has / when is"
    "when've" = "when have",
    "where'd" = "where did",
    "where's" = "where is", #"where has / where is"
    "where've" = "where have",
    "who'll" = "who will", #"who shall / who will"
    "who'll've" = "who will have", #"who shall have / who will have"
    "who's" = "who is", #"who has / who is"
    "who've" = "who have",
    "why's" = "why is", #"why has / why is"
    "why've" = "why have",
    "will've" = "will have",
    "won't" = "will not",
    "won't've" = "will not have",
    "would've" = "would have",
    "wouldn't" = "would not",
    "wouldn't've" = "would not have",
    "y'all" = "you all",
    "y'all'd" = "you all would",
    "y'all'd've" = "you all would have",
    "y'all're" = "you all are",
    "y'all've" = "you all have",
    "you'd" = "you would", #"you had / you would"
    "you'd've" = "you would have",
    "you'll" = "you will", #"you shall / you will"
    "you'll've" = "you will have", #"you shall have / you will have"
    "you're" = "you are",
    "you've" = "you have")

# data-preprocessing --------------------------------------------------------------------

# load data
data_raw <- read_csv("your_transcript.csv")
# the formatting of data should look like this
# tribble(
#   ~turn, ~content,
#   "T1", "Hello, world.",
#   "C1", "Hi, there.",
#   "T2", "This is an example.",
#   "C2", "I understand."
# )

# you can use your own list of stop words with this code. 
stopwords <- read_csv("stoplist.csv", col_names = FALSE) %>% 
  rename(word = X1)
# alternatively, if you don't have one, you can use 'stopwords <- tibble(word = stopwords::stopwords())'

# make it lowercase & replace contractions
# if you have trouble replacing contractions here, replace them in your raw data (e.g., csv file) before importing it.
data_lowercased <- data_raw %>% 
  mutate(turn = row_number()) %>% 
  mutate(content = str_to_lower(content)) %>%  #though 'unnest_tokens' make it lowercase, it needs to be done before then, to replace contractions
  mutate(content = str_replace_all(content, "¡¯", "'")) %>% 
  mutate(content = str_replace_all(content, contractions_list))

# make it into a tidytext format
data_tidied <- data_lowercased %>%
  unnest_tokens(sent, content, token = "regex", pattern = "[\\.\\?]") %>% 
  arrange(turn) %>% 
  mutate(sent_num = row_number()) %>% # give numbers for sentence window
  unnest_tokens(word, sent) %>%
  mutate(word = lemmatize_strings(word)) %>%
  filter(!str_detect(word, "[:digit:]")) %>% 
  anti_join(stopwords) %>%
  arrange(sent_num) %>% 
  mutate(sent_num = as.numeric(as.factor(sent_num))) %>% 
  mutate(win_num = (sent_num + 1) %/% 2) # sentence window size = 2, in this case // adjust it if needed

# get feature matrix
feature_matrix <- get_feature_matrix(data_tidied)

# raw turn-to-turn conceptual similarity values. Use this to calculate MPR metrics
turn_sim_mpr <- get_turn_sim(feature_matrix)


# plotting ----------------------------------------------------------------

# nonlinear transformation for visual clarity
turn_sim <- turn_sim_mpr^2 

# calculate turn_length to generate different sizes of cells. You can also use it for length filter
turn_num <- tibble(
  turn = c(1:nrow(data_raw))
)
turn_length <- data_tidied %>% 
  select(turn, word) %>% 
  distinct() %>% 
  group_by(turn) %>% 
  count() %>% 
  full_join(turn_num) %>% 
  arrange(turn) %>% 
  ungroup()
turn_length[is.na(turn_length)] <- 1
# for varying cell size in the plot
turn_length_size <- turn_length %>% 
  mutate(cs = cumsum(n)) %>% 
  mutate(posit = cs - n/2)

# change the data into long format to generate a plot
# add position and size information for plotting
turn_sim_plotting <- melt(turn_sim) %>% 
  mutate(value = ifelse(value < 0.2, 0, value)) %>% # set a cutoff value to see clearer pattern
  mutate(category = ifelse(value == 0, "N", ifelse(Var2 %% 2 == 1 & (Var1 + Var2) %% 2 == 0, "T", 
                                                   ifelse(Var2 %% 2 == 0 & (Var1 + Var2) %% 2 == 0, "C", "TC")))) %>% 
  mutate(category = factor(category)) %>% # Attach size/location information for varying size of cells 
  left_join(select(turn_length_size, turn, posit), by = c("Var2" = "turn")) %>% 
  rename(x.posit = posit) %>% 
  left_join(select(turn_length_size, turn, posit), by = c("Var1" = "turn")) %>% 
  rename(y.posit = posit) %>% 
  left_join(select(turn_length_size, turn, n), by = c("Var2" = "turn")) %>% 
  rename(w.size = n) %>% 
  left_join(select(turn_length_size, turn, n), by = c("Var1" = "turn")) %>% 
  rename(h.size = n)

# may take some time to generate the plot
get_crp(turn_sim_plotting)

# mpr ---------------------------------------------------------------------

# Draw plots for mpr (It may be hard to read, due to too frequent fluctuation)
mpr_metrics <- bind_rows(get_mpr_metrics(turn_sim, 1), get_mpr_metrics(turn_sim, 2)) %>%
  mutate(participant = factor(participant)) %>%
  group_by(participant) %>%
  mutate(turn = row_number()) %>%
  select(participant, turn, everything())
# get_mpr_plot(mpr_metrics, ITR)

# Draw plots for mpr (aggregate 3 turns into 1 for visual clarity)
mpr_metrics_3 <- mpr_metrics %>%
  mutate(turn = (((turn - 1) %/% 3) + 1) * 3) %>% 
  group_by(participant, turn) %>% 
  summarize_all(sum)
get_mpr_plot(mpr_metrics_3, ITR) # Replace ITR to the metric of your interest

# See mean values
mpr_metrics %>% group_by(participant) %>% summarize(m = mean(LTCS)) # Replace LTCS to the metric of your interest

# Interactive plotting on mpr metrics using plotly
# requires 'plotly' package
plotly::ggplotly(get_mpr_plot(mpr_metrics, ITR))


