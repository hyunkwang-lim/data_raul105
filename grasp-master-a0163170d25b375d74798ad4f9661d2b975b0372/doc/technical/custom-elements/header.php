<?php
  require($_SERVER["DOCUMENT_ROOT"] .  "/wp-load.php");
  $wp->init(); 
  $wp->parse_request(); 
  $wp->query_posts();
  $wp->register_globals(); 
  $wp->send_headers();

  class MyPost { var $post_title = "GRASP Framework Overview"; }
  $wp_query->is_home = false;
  $wp_query->is_single = true;
  $wp_query->queried_object = new MyPost();

  get_header(); ?>

<link href="$relpath^stylesheet.php.css" rel="stylesheet" type="text/css"/>
<link href="$relpath^tabs.php.css" rel="stylesheet" type="text/css"/>
<!-- <script type="text/javascript" src="$relpath^jquery.js"></script> -->
<script type="text/javascript" src="$relpath^dynsections.js"></script>
$treeview
$mathjax
<script type="text/javascript">
  window.$ = jQuery;
</script>

<div id="content" class="l-submain-h" role="main">
  <div style="padding-top: 40px;" class="l-submain">
    <div id="top"><!-- do not remove this div, it is closed by doxygen! -->