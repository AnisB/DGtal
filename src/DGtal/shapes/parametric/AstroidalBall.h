  


<!DOCTYPE html>
<html>
  <head prefix="og: http://ogp.me/ns# fb: http://ogp.me/ns/fb# githubog: http://ogp.me/ns/fb/githubog#">
    <meta charset='utf-8'>
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
        <title>DGtal/src/DGtal/shapes/parametric/AstroidalBall.h at master · AnisB/DGtal</title>
    <link rel="search" type="application/opensearchdescription+xml" href="/opensearch.xml" title="GitHub" />
    <link rel="fluid-icon" href="https://github.com/fluidicon.png" title="GitHub" />
    <link rel="shortcut icon" href="/favicon.ico" type="image/x-icon" />

    
    

    <meta content="authenticity_token" name="csrf-param" />
<meta content="yWpMBfeHy3MErLBfVmKFbhiVdg8Ff1soqJtcQb7Dzdo=" name="csrf-token" />

    <link href="https://a248.e.akamai.net/assets.github.com/stylesheets/bundles/github-83f6d9de9d84b203707ab4ba8bdf0bd6478f0d5b.css" media="screen" rel="stylesheet" type="text/css" />
    <link href="https://a248.e.akamai.net/assets.github.com/stylesheets/bundles/github2-84fb3a68186ca69e4930d889b6a55de7fb79856a.css" media="screen" rel="stylesheet" type="text/css" />
    
    


    <script src="https://a248.e.akamai.net/assets.github.com/javascripts/bundles/frameworks-a450c7f907bdc1ee6b362ab1ecca811c761fd259.js" type="text/javascript"></script>
    
    <script defer="defer" src="https://a248.e.akamai.net/assets.github.com/javascripts/bundles/github-8c27679ce3c4109571a8ad967a3cb83df487a71b.js" type="text/javascript"></script>
    
    

      <link rel='permalink' href='/AnisB/DGtal/blob/89108d7230c1510e0faf111f9128bb11a97b89e7/src/DGtal/shapes/parametric/AstroidalBall.h'>
    <meta property="og:title" content="DGtal"/>
    <meta property="og:type" content="githubog:gitrepository"/>
    <meta property="og:url" content="https://github.com/AnisB/DGtal"/>
    <meta property="og:image" content="https://a248.e.akamai.net/assets.github.com/images/gravatars/gravatar-140.png?1338942895"/>
    <meta property="og:site_name" content="GitHub"/>
    <meta property="og:description" content="Digital Geometry Tools and Algorithm Library. Contribute to DGtal development by creating an account on GitHub."/>

    <meta name="description" content="Digital Geometry Tools and Algorithm Library. Contribute to DGtal development by creating an account on GitHub." />

  <link href="https://github.com/AnisB/DGtal/commits/master.atom" rel="alternate" title="Recent Commits to DGtal:master" type="application/atom+xml" />

  </head>


  <body class="logged_in page-blob linux vis-public fork env-production " data-blob-contribs-enabled="yes">
    <div id="wrapper">

    
    
    

      <div id="header" class="true clearfix">
        <div class="container clearfix">
          <a class="site-logo" href="https://github.com/">
            <!--[if IE]>
            <img alt="GitHub" class="github-logo" src="https://a248.e.akamai.net/assets.github.com/images/modules/header/logov7.png?1338942895" />
            <img alt="GitHub" class="github-logo-hover" src="https://a248.e.akamai.net/assets.github.com/images/modules/header/logov7-hover.png?1338942895" />
            <![endif]-->
            <img alt="GitHub" class="github-logo-4x" height="30" src="https://a248.e.akamai.net/assets.github.com/images/modules/header/logov7@4x.png?1338942895" />
            <img alt="GitHub" class="github-logo-4x-hover" height="30" src="https://a248.e.akamai.net/assets.github.com/images/modules/header/logov7@4x-hover.png?1338942895" />
          </a>

              
    <div class="topsearch  ">
        <form accept-charset="UTF-8" action="/search" id="top_search_form" method="get">
  <a href="/search" class="advanced-search tooltipped downwards" title="Advanced Search"><span class="mini-icon mini-icon-advanced-search"></span></a>
  <div class="search placeholder-field js-placeholder-field">
    <input type="text" class="search my_repos_autocompleter" id="global-search-field" name="q" results="5" spellcheck="false" autocomplete="off" data-autocomplete="my-repos-autocomplete" placeholder="Search&hellip;" data-hotkey="s" />
    <div id="my-repos-autocomplete" class="autocomplete-results">
      <ul class="js-navigation-container"></ul>
    </div>
    <input type="submit" value="Search" class="button">
    <span class="mini-icon mini-icon-search-input"></span>
  </div>
  <input type="hidden" name="type" value="Everything" />
  <input type="hidden" name="repo" value="" />
  <input type="hidden" name="langOverride" value="" />
  <input type="hidden" name="start_value" value="1" />
</form>

      <ul class="top-nav">
          <li class="explore"><a href="https://github.com/explore">Explore</a></li>
          <li><a href="https://gist.github.com">Gist</a></li>
          <li><a href="/blog">Blog</a></li>
        <li><a href="http://help.github.com">Help</a></li>
      </ul>
    </div>


            


  <div id="userbox">
    <div id="user">
      <a href="https://github.com/AnisB"><img height="20" src="https://secure.gravatar.com/avatar/3d1b91e4a19046ddf1ce905715acd71d?s=140&amp;d=https://a248.e.akamai.net/assets.github.com%2Fimages%2Fgravatars%2Fgravatar-140.png" width="20" /></a>
      <a href="https://github.com/AnisB" class="name">AnisB</a>
    </div>
    <ul id="user-links">
      <li>
        <a href="/new" id="new_repo" class="tooltipped downwards" title="Create a New Repo"><span class="mini-icon mini-icon-create"></span></a>
      </li>
      <li>
        <a href="/inbox/notifications" id="notifications" class="tooltipped downwards" title="Notifications">
          <span class="mini-icon mini-icon-notifications"></span>
        </a>
      </li>
      <li>
        <a href="/settings/profile" id="settings" class="tooltipped downwards" title="Account Settings">
          <span class="mini-icon mini-icon-account-settings"></span>
        </a>
      </li>
      <li>
          <a href="/logout" data-method="post" id="logout" class="tooltipped downwards" title="Log Out">
            <span class="mini-icon mini-icon-logout"></span>
          </a>
      </li>
    </ul>
  </div>



          
        </div>
      </div>

      

            <div class="site hfeed" itemscope itemtype="http://schema.org/WebPage">
      <div class="container hentry">
        <div class="pagehead repohead instapaper_ignore readability-menu">
        <div class="title-actions-bar">
          



              <ul class="pagehead-actions">

          <li class="for-owner"><a href="/AnisB/DGtal/admin" class="minibutton btn-admin "><span class="icon"></span>Admin</a></li>

              <li class="nspr">
                <a href="/AnisB/DGtal/pull/new/master" class="minibutton btn-pull-request "><span class="icon"></span>Pull Request</a>
              </li>

          <li class="js-toggler-container js-social-container watch-button-container ">
            <span class="watch-button"><a href="/AnisB/DGtal/toggle_watch" class="minibutton btn-watch js-toggler-target lighter" data-remote="true" data-method="post" rel="nofollow"><span class="icon"></span> Watch</a><a class="social-count js-social-count" href="/AnisB/DGtal/watchers">0</a></span>
            <span class="unwatch-button"><a href="/AnisB/DGtal/toggle_watch" class="minibutton btn-unwatch js-toggler-target lighter" data-remote="true" data-method="post" rel="nofollow"><span class="icon"></span> Unwatch</a><a class="social-count js-social-count" href="/AnisB/DGtal/watchers">0</a></span>
          </li>

              <li><a href="/AnisB/DGtal/fork" class="minibutton btn-fork js-toggler-target fork-button lighter" rel="nofollow" data-method="post"><span class="icon"></span> Fork</a><a href="/AnisB/DGtal/network" class="social-count">18</a>
              </li>


    </ul>

          <h1 itemscope itemtype="http://data-vocabulary.org/Breadcrumb" class="entry-title public">
            <span class="repo-label"><span>public</span></span>
            <span class="mega-icon mega-icon-repo-forked"></span>
            <span class="author vcard">
<a href="/AnisB" class="url fn" itemprop="url" rel="author">              <span itemprop="title">AnisB</span>
              </a></span> /
            <strong><a href="/AnisB/DGtal" class="js-current-repository">DGtal</a></strong>
              <span class="fork-flag">
                <span class="text">forked from <a href="/DGtal-team/DGtal">DGtal-team/DGtal</a></span>
              </span>
          </h1>
        </div>

          

  <ul class="tabs">
    <li><a href="/AnisB/DGtal" class="selected" highlight="repo_sourcerepo_downloadsrepo_commitsrepo_tagsrepo_branches">Code</a></li>
    <li><a href="/AnisB/DGtal/network" highlight="repo_network">Network</a>
    <li><a href="/AnisB/DGtal/pulls" highlight="repo_pulls">Pull Requests <span class='counter'>0</span></a></li>


      <li><a href="/AnisB/DGtal/wiki" highlight="repo_wiki">Wiki</a></li>

    <li><a href="/AnisB/DGtal/graphs" highlight="repo_graphsrepo_contributors">Graphs</a></li>

  </ul>
 
<div class="frame frame-center tree-finder" style="display:none"
      data-tree-list-url="/AnisB/DGtal/tree-list/89108d7230c1510e0faf111f9128bb11a97b89e7"
      data-blob-url-prefix="/AnisB/DGtal/blob/89108d7230c1510e0faf111f9128bb11a97b89e7"
    >

  <div class="breadcrumb">
    <span class="bold"><a href="/AnisB/DGtal">DGtal</a></span> /
    <input class="tree-finder-input js-navigation-enable" type="text" name="query" autocomplete="off" spellcheck="false">
  </div>

    <div class="octotip">
      <p>
        <a href="/AnisB/DGtal/dismiss-tree-finder-help" class="dismiss js-dismiss-tree-list-help" title="Hide this notice forever" rel="nofollow">Dismiss</a>
        <span class="bold">Octotip:</span> You've activated the <em>file finder</em>
        by pressing <span class="kbd">t</span> Start typing to filter the
        file list. Use <span class="kbd badmono">↑</span> and
        <span class="kbd badmono">↓</span> to navigate,
        <span class="kbd">enter</span> to view files.
      </p>
    </div>

  <table class="tree-browser" cellpadding="0" cellspacing="0">
    <tr class="js-header"><th>&nbsp;</th><th>name</th></tr>
    <tr class="js-no-results no-results" style="display: none">
      <th colspan="2">No matching files</th>
    </tr>
    <tbody class="js-results-list js-navigation-container">
    </tbody>
  </table>
</div>

<div id="jump-to-line" style="display:none">
  <h2>Jump to Line</h2>
  <form accept-charset="UTF-8">
    <input class="textfield" type="text">
    <div class="full-button">
      <button type="submit" class="classy">
        <span>Go</span>
      </button>
    </div>
  </form>
</div>


<div class="subnav-bar">

  <ul class="actions subnav">
    <li><a href="/AnisB/DGtal/tags" class="" highlight="repo_tags">Tags <span class="counter">3</span></a></li>
    <li><a href="/AnisB/DGtal/downloads" class="blank downloads-blank" highlight="repo_downloads">Downloads <span class="counter">0</span></a></li>
    
  </ul>

  <ul class="scope">
    <li class="switcher">

      <div class="context-menu-container js-menu-container js-context-menu">
        <a href="#"
           class="minibutton bigger switcher js-menu-target js-commitish-button btn-branch repo-tree"
           data-hotkey="w"
           data-master-branch="master"
           data-ref="master">
           <span><span class="icon"></span><i>branch:</i> master</span>
        </a>

        <div class="context-pane commitish-context js-menu-content">
          <a href="javascript:;" class="close js-menu-close"><span class="mini-icon mini-icon-remove-close"></span></a>
          <div class="context-title">Switch Branches/Tags</div>
          <div class="context-body pane-selector commitish-selector js-navigation-container">
            <div class="filterbar">
              <input type="text" id="context-commitish-filter-field" class="js-navigation-enable" placeholder="Filter branches/tags" data-filterable />

              <ul class="tabs">
                <li><a href="#" data-filter="branches" class="selected">Branches</a></li>
                <li><a href="#" data-filter="tags">Tags</a></li>
              </ul>
            </div>

            <div class="js-filter-tab js-filter-branches" data-filterable-for="context-commitish-filter-field">
              <div class="no-results js-not-filterable">Nothing to show</div>
                <div class="commitish-item branch-commitish selector-item js-navigation-item js-navigation-target">
                  <h4>
                      <a href="/AnisB/DGtal/blob/dcoeurjo/master/src/DGtal/shapes/parametric/AstroidalBall.h" class="js-navigation-open" data-name="dcoeurjo/master" rel="nofollow">dcoeurjo/master</a>
                  </h4>
                </div>
                <div class="commitish-item branch-commitish selector-item js-navigation-item js-navigation-target">
                  <h4>
                      <a href="/AnisB/DGtal/blob/digitalSnow/src/DGtal/shapes/parametric/AstroidalBall.h" class="js-navigation-open" data-name="digitalSnow" rel="nofollow">digitalSnow</a>
                  </h4>
                </div>
                <div class="commitish-item branch-commitish selector-item js-navigation-item js-navigation-target">
                  <h4>
                      <a href="/AnisB/DGtal/blob/master/src/DGtal/shapes/parametric/AstroidalBall.h" class="js-navigation-open" data-name="master" rel="nofollow">master</a>
                  </h4>
                </div>
                <div class="commitish-item branch-commitish selector-item js-navigation-item js-navigation-target">
                  <h4>
                      <a href="/AnisB/DGtal/blob/ViewerParam/src/DGtal/shapes/parametric/AstroidalBall.h" class="js-navigation-open" data-name="ViewerParam" rel="nofollow">ViewerParam</a>
                  </h4>
                </div>
            </div>

            <div class="js-filter-tab js-filter-tags" style="display:none" data-filterable-for="context-commitish-filter-field">
              <div class="no-results js-not-filterable">Nothing to show</div>
                <div class="commitish-item tag-commitish selector-item js-navigation-item js-navigation-target">
                  <h4>
                      <a href="/AnisB/DGtal/blob/v0.5.1/src/DGtal/shapes/parametric/AstroidalBall.h" class="js-navigation-open" data-name="v0.5.1" rel="nofollow">v0.5.1</a>
                  </h4>
                </div>
                <div class="commitish-item tag-commitish selector-item js-navigation-item js-navigation-target">
                  <h4>
                      <a href="/AnisB/DGtal/blob/v0.5/src/DGtal/shapes/parametric/AstroidalBall.h" class="js-navigation-open" data-name="v0.5" rel="nofollow">v0.5</a>
                  </h4>
                </div>
                <div class="commitish-item tag-commitish selector-item js-navigation-item js-navigation-target">
                  <h4>
                      <a href="/AnisB/DGtal/blob/0.4/src/DGtal/shapes/parametric/AstroidalBall.h" class="js-navigation-open" data-name="0.4" rel="nofollow">0.4</a>
                  </h4>
                </div>
            </div>
          </div>
        </div><!-- /.commitish-context-context -->
      </div>

    </li>
  </ul>

  <ul class="subnav with-scope">

    <li><a href="/AnisB/DGtal" class="selected" highlight="repo_source">Files</a></li>
    <li><a href="/AnisB/DGtal/commits/master" highlight="repo_commits">Commits</a></li>
    <li><a href="/AnisB/DGtal/branches" class="" highlight="repo_branches" rel="nofollow">Branches <span class="counter">4</span></a></li>
  </ul>

</div>

  
  
  


          

        </div><!-- /.repohead -->

        





<!-- block_view_fragment_key: views10/v8/blob:v21:2dba9b67ad1fc6102ace788f7281f3bb -->
  <div id="slider">

    <div class="breadcrumb" data-path="src/DGtal/shapes/parametric/AstroidalBall.h/">
      <b itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/AnisB/DGtal/tree/89108d7230c1510e0faf111f9128bb11a97b89e7" class="js-rewrite-sha" itemprop="url"><span itemprop="title">DGtal</span></a></b> / <span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/AnisB/DGtal/tree/89108d7230c1510e0faf111f9128bb11a97b89e7/src" class="js-rewrite-sha" itemscope="url"><span itemprop="title">src</span></a></span> / <span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/AnisB/DGtal/tree/89108d7230c1510e0faf111f9128bb11a97b89e7/src/DGtal" class="js-rewrite-sha" itemscope="url"><span itemprop="title">DGtal</span></a></span> / <span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/AnisB/DGtal/tree/89108d7230c1510e0faf111f9128bb11a97b89e7/src/DGtal/shapes" class="js-rewrite-sha" itemscope="url"><span itemprop="title">shapes</span></a></span> / <span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/AnisB/DGtal/tree/89108d7230c1510e0faf111f9128bb11a97b89e7/src/DGtal/shapes/parametric" class="js-rewrite-sha" itemscope="url"><span itemprop="title">parametric</span></a></span> / <strong class="final-path">AstroidalBall.h</strong> <span class="js-clippy mini-icon mini-icon-clippy " data-clipboard-text="src/DGtal/shapes/parametric/AstroidalBall.h" data-copied-hint="copied!" data-copy-hint="copy to clipboard"></span>
    </div>


      <div class="commit file-history-tease" data-path="src/DGtal/shapes/parametric/AstroidalBall.h/">
        <img class="main-avatar" height="24" src="https://secure.gravatar.com/avatar/25b38440b41b8da4492117a828121d79?s=140&amp;d=https://a248.e.akamai.net/assets.github.com%2Fimages%2Fgravatars%2Fgravatar-140.png" width="24" />
        <span class="author">Anis</span>
        <time class="js-relative-date" datetime="2012-06-11T03:55:33-07:00" title="2012-06-11 03:55:33">June 11, 2012</time>
        <div class="commit-title">
            <a href="/AnisB/DGtal/commit/f998b4654cdf60e5d853cb5a74c7627b1c701a1d" class="message">Commit starshaped3D</a>
        </div>

        <div class="participation">
          <p class="quickstat"><a href="#blob_contributors_box" rel="facebox"><strong>0</strong> contributors</a></p>
          
        </div>
        <div id="blob_contributors_box" style="display:none">
          <h2>Users on GitHub who have contributed to this file</h2>
          <ul class="facebox-user-list">
          </ul>
        </div>
      </div>

    <div class="frames">
      <div class="frame frame-center" data-path="src/DGtal/shapes/parametric/AstroidalBall.h/" data-permalink-url="/AnisB/DGtal/blob/89108d7230c1510e0faf111f9128bb11a97b89e7/src/DGtal/shapes/parametric/AstroidalBall.h" data-title="DGtal/src/DGtal/shapes/parametric/AstroidalBall.h at master · AnisB/DGtal · GitHub" data-type="blob">

        <div id="files" class="bubble">
          <div class="file">
            <div class="meta">
              <div class="info">
                <span class="icon"><b class="mini-icon mini-icon-text-file"></b></span>
                <span class="mode" title="File Mode">100644</span>
                  <span>292 lines (225 sloc)</span>
                <span>8.272 kb</span>
              </div>
              <ul class="button-group actions">
                  <li>
                    <a class="grouped-button file-edit-link minibutton bigger lighter js-rewrite-sha" href="/AnisB/DGtal/edit/89108d7230c1510e0faf111f9128bb11a97b89e7/src/DGtal/shapes/parametric/AstroidalBall.h" data-method="post" rel="nofollow" data-hotkey="e">Edit this file</a>
                  </li>

                <li>
                  <a href="/AnisB/DGtal/raw/master/src/DGtal/shapes/parametric/AstroidalBall.h" class="minibutton btn-raw grouped-button bigger lighter" id="raw-url"><span class="icon"></span>Raw</a>
                </li>
                  <li>
                    <a href="/AnisB/DGtal/blame/master/src/DGtal/shapes/parametric/AstroidalBall.h" class="minibutton btn-blame grouped-button bigger lighter"><span class="icon"></span>Blame</a>
                  </li>
                <li>
                  <a href="/AnisB/DGtal/commits/master/src/DGtal/shapes/parametric/AstroidalBall.h" class="minibutton btn-history grouped-button bigger lighter" rel="nofollow"><span class="icon"></span>History</a>
                </li>
              </ul>
            </div>
              <div class="data type-c">
      <table cellpadding="0" cellspacing="0" class="lines">
        <tr>
          <td>
            <pre class="line_numbers"><span id="L1" rel="#L1">1</span>
<span id="L2" rel="#L2">2</span>
<span id="L3" rel="#L3">3</span>
<span id="L4" rel="#L4">4</span>
<span id="L5" rel="#L5">5</span>
<span id="L6" rel="#L6">6</span>
<span id="L7" rel="#L7">7</span>
<span id="L8" rel="#L8">8</span>
<span id="L9" rel="#L9">9</span>
<span id="L10" rel="#L10">10</span>
<span id="L11" rel="#L11">11</span>
<span id="L12" rel="#L12">12</span>
<span id="L13" rel="#L13">13</span>
<span id="L14" rel="#L14">14</span>
<span id="L15" rel="#L15">15</span>
<span id="L16" rel="#L16">16</span>
<span id="L17" rel="#L17">17</span>
<span id="L18" rel="#L18">18</span>
<span id="L19" rel="#L19">19</span>
<span id="L20" rel="#L20">20</span>
<span id="L21" rel="#L21">21</span>
<span id="L22" rel="#L22">22</span>
<span id="L23" rel="#L23">23</span>
<span id="L24" rel="#L24">24</span>
<span id="L25" rel="#L25">25</span>
<span id="L26" rel="#L26">26</span>
<span id="L27" rel="#L27">27</span>
<span id="L28" rel="#L28">28</span>
<span id="L29" rel="#L29">29</span>
<span id="L30" rel="#L30">30</span>
<span id="L31" rel="#L31">31</span>
<span id="L32" rel="#L32">32</span>
<span id="L33" rel="#L33">33</span>
<span id="L34" rel="#L34">34</span>
<span id="L35" rel="#L35">35</span>
<span id="L36" rel="#L36">36</span>
<span id="L37" rel="#L37">37</span>
<span id="L38" rel="#L38">38</span>
<span id="L39" rel="#L39">39</span>
<span id="L40" rel="#L40">40</span>
<span id="L41" rel="#L41">41</span>
<span id="L42" rel="#L42">42</span>
<span id="L43" rel="#L43">43</span>
<span id="L44" rel="#L44">44</span>
<span id="L45" rel="#L45">45</span>
<span id="L46" rel="#L46">46</span>
<span id="L47" rel="#L47">47</span>
<span id="L48" rel="#L48">48</span>
<span id="L49" rel="#L49">49</span>
<span id="L50" rel="#L50">50</span>
<span id="L51" rel="#L51">51</span>
<span id="L52" rel="#L52">52</span>
<span id="L53" rel="#L53">53</span>
<span id="L54" rel="#L54">54</span>
<span id="L55" rel="#L55">55</span>
<span id="L56" rel="#L56">56</span>
<span id="L57" rel="#L57">57</span>
<span id="L58" rel="#L58">58</span>
<span id="L59" rel="#L59">59</span>
<span id="L60" rel="#L60">60</span>
<span id="L61" rel="#L61">61</span>
<span id="L62" rel="#L62">62</span>
<span id="L63" rel="#L63">63</span>
<span id="L64" rel="#L64">64</span>
<span id="L65" rel="#L65">65</span>
<span id="L66" rel="#L66">66</span>
<span id="L67" rel="#L67">67</span>
<span id="L68" rel="#L68">68</span>
<span id="L69" rel="#L69">69</span>
<span id="L70" rel="#L70">70</span>
<span id="L71" rel="#L71">71</span>
<span id="L72" rel="#L72">72</span>
<span id="L73" rel="#L73">73</span>
<span id="L74" rel="#L74">74</span>
<span id="L75" rel="#L75">75</span>
<span id="L76" rel="#L76">76</span>
<span id="L77" rel="#L77">77</span>
<span id="L78" rel="#L78">78</span>
<span id="L79" rel="#L79">79</span>
<span id="L80" rel="#L80">80</span>
<span id="L81" rel="#L81">81</span>
<span id="L82" rel="#L82">82</span>
<span id="L83" rel="#L83">83</span>
<span id="L84" rel="#L84">84</span>
<span id="L85" rel="#L85">85</span>
<span id="L86" rel="#L86">86</span>
<span id="L87" rel="#L87">87</span>
<span id="L88" rel="#L88">88</span>
<span id="L89" rel="#L89">89</span>
<span id="L90" rel="#L90">90</span>
<span id="L91" rel="#L91">91</span>
<span id="L92" rel="#L92">92</span>
<span id="L93" rel="#L93">93</span>
<span id="L94" rel="#L94">94</span>
<span id="L95" rel="#L95">95</span>
<span id="L96" rel="#L96">96</span>
<span id="L97" rel="#L97">97</span>
<span id="L98" rel="#L98">98</span>
<span id="L99" rel="#L99">99</span>
<span id="L100" rel="#L100">100</span>
<span id="L101" rel="#L101">101</span>
<span id="L102" rel="#L102">102</span>
<span id="L103" rel="#L103">103</span>
<span id="L104" rel="#L104">104</span>
<span id="L105" rel="#L105">105</span>
<span id="L106" rel="#L106">106</span>
<span id="L107" rel="#L107">107</span>
<span id="L108" rel="#L108">108</span>
<span id="L109" rel="#L109">109</span>
<span id="L110" rel="#L110">110</span>
<span id="L111" rel="#L111">111</span>
<span id="L112" rel="#L112">112</span>
<span id="L113" rel="#L113">113</span>
<span id="L114" rel="#L114">114</span>
<span id="L115" rel="#L115">115</span>
<span id="L116" rel="#L116">116</span>
<span id="L117" rel="#L117">117</span>
<span id="L118" rel="#L118">118</span>
<span id="L119" rel="#L119">119</span>
<span id="L120" rel="#L120">120</span>
<span id="L121" rel="#L121">121</span>
<span id="L122" rel="#L122">122</span>
<span id="L123" rel="#L123">123</span>
<span id="L124" rel="#L124">124</span>
<span id="L125" rel="#L125">125</span>
<span id="L126" rel="#L126">126</span>
<span id="L127" rel="#L127">127</span>
<span id="L128" rel="#L128">128</span>
<span id="L129" rel="#L129">129</span>
<span id="L130" rel="#L130">130</span>
<span id="L131" rel="#L131">131</span>
<span id="L132" rel="#L132">132</span>
<span id="L133" rel="#L133">133</span>
<span id="L134" rel="#L134">134</span>
<span id="L135" rel="#L135">135</span>
<span id="L136" rel="#L136">136</span>
<span id="L137" rel="#L137">137</span>
<span id="L138" rel="#L138">138</span>
<span id="L139" rel="#L139">139</span>
<span id="L140" rel="#L140">140</span>
<span id="L141" rel="#L141">141</span>
<span id="L142" rel="#L142">142</span>
<span id="L143" rel="#L143">143</span>
<span id="L144" rel="#L144">144</span>
<span id="L145" rel="#L145">145</span>
<span id="L146" rel="#L146">146</span>
<span id="L147" rel="#L147">147</span>
<span id="L148" rel="#L148">148</span>
<span id="L149" rel="#L149">149</span>
<span id="L150" rel="#L150">150</span>
<span id="L151" rel="#L151">151</span>
<span id="L152" rel="#L152">152</span>
<span id="L153" rel="#L153">153</span>
<span id="L154" rel="#L154">154</span>
<span id="L155" rel="#L155">155</span>
<span id="L156" rel="#L156">156</span>
<span id="L157" rel="#L157">157</span>
<span id="L158" rel="#L158">158</span>
<span id="L159" rel="#L159">159</span>
<span id="L160" rel="#L160">160</span>
<span id="L161" rel="#L161">161</span>
<span id="L162" rel="#L162">162</span>
<span id="L163" rel="#L163">163</span>
<span id="L164" rel="#L164">164</span>
<span id="L165" rel="#L165">165</span>
<span id="L166" rel="#L166">166</span>
<span id="L167" rel="#L167">167</span>
<span id="L168" rel="#L168">168</span>
<span id="L169" rel="#L169">169</span>
<span id="L170" rel="#L170">170</span>
<span id="L171" rel="#L171">171</span>
<span id="L172" rel="#L172">172</span>
<span id="L173" rel="#L173">173</span>
<span id="L174" rel="#L174">174</span>
<span id="L175" rel="#L175">175</span>
<span id="L176" rel="#L176">176</span>
<span id="L177" rel="#L177">177</span>
<span id="L178" rel="#L178">178</span>
<span id="L179" rel="#L179">179</span>
<span id="L180" rel="#L180">180</span>
<span id="L181" rel="#L181">181</span>
<span id="L182" rel="#L182">182</span>
<span id="L183" rel="#L183">183</span>
<span id="L184" rel="#L184">184</span>
<span id="L185" rel="#L185">185</span>
<span id="L186" rel="#L186">186</span>
<span id="L187" rel="#L187">187</span>
<span id="L188" rel="#L188">188</span>
<span id="L189" rel="#L189">189</span>
<span id="L190" rel="#L190">190</span>
<span id="L191" rel="#L191">191</span>
<span id="L192" rel="#L192">192</span>
<span id="L193" rel="#L193">193</span>
<span id="L194" rel="#L194">194</span>
<span id="L195" rel="#L195">195</span>
<span id="L196" rel="#L196">196</span>
<span id="L197" rel="#L197">197</span>
<span id="L198" rel="#L198">198</span>
<span id="L199" rel="#L199">199</span>
<span id="L200" rel="#L200">200</span>
<span id="L201" rel="#L201">201</span>
<span id="L202" rel="#L202">202</span>
<span id="L203" rel="#L203">203</span>
<span id="L204" rel="#L204">204</span>
<span id="L205" rel="#L205">205</span>
<span id="L206" rel="#L206">206</span>
<span id="L207" rel="#L207">207</span>
<span id="L208" rel="#L208">208</span>
<span id="L209" rel="#L209">209</span>
<span id="L210" rel="#L210">210</span>
<span id="L211" rel="#L211">211</span>
<span id="L212" rel="#L212">212</span>
<span id="L213" rel="#L213">213</span>
<span id="L214" rel="#L214">214</span>
<span id="L215" rel="#L215">215</span>
<span id="L216" rel="#L216">216</span>
<span id="L217" rel="#L217">217</span>
<span id="L218" rel="#L218">218</span>
<span id="L219" rel="#L219">219</span>
<span id="L220" rel="#L220">220</span>
<span id="L221" rel="#L221">221</span>
<span id="L222" rel="#L222">222</span>
<span id="L223" rel="#L223">223</span>
<span id="L224" rel="#L224">224</span>
<span id="L225" rel="#L225">225</span>
<span id="L226" rel="#L226">226</span>
<span id="L227" rel="#L227">227</span>
<span id="L228" rel="#L228">228</span>
<span id="L229" rel="#L229">229</span>
<span id="L230" rel="#L230">230</span>
<span id="L231" rel="#L231">231</span>
<span id="L232" rel="#L232">232</span>
<span id="L233" rel="#L233">233</span>
<span id="L234" rel="#L234">234</span>
<span id="L235" rel="#L235">235</span>
<span id="L236" rel="#L236">236</span>
<span id="L237" rel="#L237">237</span>
<span id="L238" rel="#L238">238</span>
<span id="L239" rel="#L239">239</span>
<span id="L240" rel="#L240">240</span>
<span id="L241" rel="#L241">241</span>
<span id="L242" rel="#L242">242</span>
<span id="L243" rel="#L243">243</span>
<span id="L244" rel="#L244">244</span>
<span id="L245" rel="#L245">245</span>
<span id="L246" rel="#L246">246</span>
<span id="L247" rel="#L247">247</span>
<span id="L248" rel="#L248">248</span>
<span id="L249" rel="#L249">249</span>
<span id="L250" rel="#L250">250</span>
<span id="L251" rel="#L251">251</span>
<span id="L252" rel="#L252">252</span>
<span id="L253" rel="#L253">253</span>
<span id="L254" rel="#L254">254</span>
<span id="L255" rel="#L255">255</span>
<span id="L256" rel="#L256">256</span>
<span id="L257" rel="#L257">257</span>
<span id="L258" rel="#L258">258</span>
<span id="L259" rel="#L259">259</span>
<span id="L260" rel="#L260">260</span>
<span id="L261" rel="#L261">261</span>
<span id="L262" rel="#L262">262</span>
<span id="L263" rel="#L263">263</span>
<span id="L264" rel="#L264">264</span>
<span id="L265" rel="#L265">265</span>
<span id="L266" rel="#L266">266</span>
<span id="L267" rel="#L267">267</span>
<span id="L268" rel="#L268">268</span>
<span id="L269" rel="#L269">269</span>
<span id="L270" rel="#L270">270</span>
<span id="L271" rel="#L271">271</span>
<span id="L272" rel="#L272">272</span>
<span id="L273" rel="#L273">273</span>
<span id="L274" rel="#L274">274</span>
<span id="L275" rel="#L275">275</span>
<span id="L276" rel="#L276">276</span>
<span id="L277" rel="#L277">277</span>
<span id="L278" rel="#L278">278</span>
<span id="L279" rel="#L279">279</span>
<span id="L280" rel="#L280">280</span>
<span id="L281" rel="#L281">281</span>
<span id="L282" rel="#L282">282</span>
<span id="L283" rel="#L283">283</span>
<span id="L284" rel="#L284">284</span>
<span id="L285" rel="#L285">285</span>
<span id="L286" rel="#L286">286</span>
<span id="L287" rel="#L287">287</span>
<span id="L288" rel="#L288">288</span>
<span id="L289" rel="#L289">289</span>
<span id="L290" rel="#L290">290</span>
<span id="L291" rel="#L291">291</span>
</pre>
          </td>
          <td width="100%">
                <div class="highlight"><pre><div class='line' id='LC1'><span class="cm">/**</span></div><div class='line' id='LC2'><span class="cm"> *  This program is free software: you can redistribute it and/or modify</span></div><div class='line' id='LC3'><span class="cm"> *  it under the terms of the GNU Lesser General Public License as</span></div><div class='line' id='LC4'><span class="cm"> *  published by the Free Software Foundation, either version 3 of the</span></div><div class='line' id='LC5'><span class="cm"> *  License, or  (at your option) any later version.</span></div><div class='line' id='LC6'><span class="cm"> *</span></div><div class='line' id='LC7'><span class="cm"> *  This program is distributed in the hope that it will be useful,</span></div><div class='line' id='LC8'><span class="cm"> *  but WITHOUT ANY WARRANTY; without even the implied warranty of</span></div><div class='line' id='LC9'><span class="cm"> *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the</span></div><div class='line' id='LC10'><span class="cm"> *  GNU General Public License for more details.</span></div><div class='line' id='LC11'><span class="cm"> *</span></div><div class='line' id='LC12'><span class="cm"> *  You should have received a copy of the GNU General Public License</span></div><div class='line' id='LC13'><span class="cm"> *  along with this program.  If not, see &lt;http://www.gnu.org/licenses/&gt;.</span></div><div class='line' id='LC14'><span class="cm"> *</span></div><div class='line' id='LC15'><span class="cm"> **/</span></div><div class='line' id='LC16'><br/></div><div class='line' id='LC17'><span class="cp">#pragma once</span></div><div class='line' id='LC18'><br/></div><div class='line' id='LC19'><span class="cm">/**</span></div><div class='line' id='LC20'><span class="cm"> * @file AstroidalBall.h</span></div><div class='line' id='LC21'><span class="cm"> * @author Anis Benyoub (\c anis.benyoub@liris.cnrs.fr )</span></div><div class='line' id='LC22'><span class="cm"> * Laboratoire d&#39;InfoRmatique en Image et Systèmes d&#39;information - LIRIS (CNRS, UMR 5205), CNRS, France</span></div><div class='line' id='LC23'><span class="cm"> *</span></div><div class='line' id='LC24'><span class="cm"> * @date 2012/06/05</span></div><div class='line' id='LC25'><span class="cm"> *</span></div><div class='line' id='LC26'><span class="cm"> * Header file for module AstroidalBall.cpp</span></div><div class='line' id='LC27'><span class="cm"> *</span></div><div class='line' id='LC28'><span class="cm"> * This file is part of the DGtal library.</span></div><div class='line' id='LC29'><span class="cm"> */</span></div><div class='line' id='LC30'><br/></div><div class='line' id='LC31'><span class="cp">#if defined(AstroidalBall_RECURSES)</span></div><div class='line' id='LC32'><span class="cp">#error Recursive header files inclusion detected in AstroidalBall.h</span></div><div class='line' id='LC33'><span class="cp">#else </span><span class="c1">// defined(AstroidalBall_RECURSES)</span></div><div class='line' id='LC34'><span class="cm">/** Prevents recursive inclusion of headers. */</span></div><div class='line' id='LC35'><span class="cp">#define AstroidalBall_RECURSES</span></div><div class='line' id='LC36'><br/></div><div class='line' id='LC37'><span class="cp">#if !defined AstroidalBall_h</span></div><div class='line' id='LC38'><span class="cm">/** Prevents repeated inclusion of headers. */</span></div><div class='line' id='LC39'><span class="cp">#define AstroidalBall_h</span></div><div class='line' id='LC40'><br/></div><div class='line' id='LC41'><span class="c1">//////////////////////////////////////////////////////////////////////////////</span></div><div class='line' id='LC42'><span class="c1">// Inclusions</span></div><div class='line' id='LC43'><span class="cp">#include &quot;DGtal/base/Common.h&quot;</span></div><div class='line' id='LC44'><span class="cp">#include &quot;DGtal/shapes/parametric/StarShaped3D.h&quot;</span></div><div class='line' id='LC45'><span class="c1">//////////////////////////////////////////////////////////////////////////////</span></div><div class='line' id='LC46'><br/></div><div class='line' id='LC47'><span class="k">namespace</span> <span class="n">DGtal</span></div><div class='line' id='LC48'><span class="p">{</span></div><div class='line' id='LC49'><br/></div><div class='line' id='LC50'>&nbsp;&nbsp;<span class="c1">/////////////////////////////////////////////////////////////////////////////</span></div><div class='line' id='LC51'>&nbsp;&nbsp;<span class="c1">// template class AstroidalBall</span></div><div class='line' id='LC52'>&nbsp;&nbsp;<span class="cm">/**</span></div><div class='line' id='LC53'><span class="cm">   * Description of template class &#39;AstroidalBall&#39; &lt;p&gt;</span></div><div class='line' id='LC54'><span class="cm">   * \brief Aim: Model of the concept StarShaped</span></div><div class='line' id='LC55'><span class="cm">   * represents any astroidalball in the space.</span></div><div class='line' id='LC56'><span class="cm">   *</span></div><div class='line' id='LC57'><span class="cm">   */</span></div><div class='line' id='LC58'>&nbsp;&nbsp;<span class="k">template</span> <span class="o">&lt;</span><span class="k">typename</span> <span class="n">TSpace</span><span class="o">&gt;</span></div><div class='line' id='LC59'>&nbsp;&nbsp;<span class="k">class</span> <span class="nc">AstroidalBall</span><span class="o">:</span>  <span class="k">public</span> <span class="n">StarShaped3D</span><span class="o">&lt;</span><span class="n">TSpace</span><span class="o">&gt;</span></div><div class='line' id='LC60'>&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC61'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// ----------------------- Standard services ------------------------------</span></div><div class='line' id='LC62'>&nbsp;&nbsp;<span class="k">public</span><span class="o">:</span></div><div class='line' id='LC63'><br/></div><div class='line' id='LC64'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">typedef</span> <span class="n">TSpace</span> <span class="n">Space</span><span class="p">;</span></div><div class='line' id='LC65'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">typedef</span> <span class="k">typename</span> <span class="n">Space</span><span class="o">::</span><span class="n">RealPoint</span> <span class="n">RealPoint</span><span class="p">;</span></div><div class='line' id='LC66'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">typedef</span>  <span class="n">pair</span><span class="o">&lt;</span><span class="kt">double</span><span class="p">,</span><span class="kt">double</span><span class="o">&gt;</span> <span class="n">AngularCoordinates</span><span class="p">;</span></div><div class='line' id='LC67'>&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC68'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="cm">/**</span></div><div class='line' id='LC69'><span class="cm">     * Destructor.</span></div><div class='line' id='LC70'><span class="cm">     */</span></div><div class='line' id='LC71'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="o">~</span><span class="n">AstroidalBall</span><span class="p">();</span></div><div class='line' id='LC72'>&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC73'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="cm">/**</span></div><div class='line' id='LC74'><span class="cm">     * Constructor. </span></div><div class='line' id='LC75'><span class="cm">     * @param x0 the x-coordinate of the AstroidalBall center.</span></div><div class='line' id='LC76'><span class="cm">     * @param y0 the y-coordinate of the AstroidalBall center.</span></div><div class='line' id='LC77'><span class="cm">     * @param y0 the y-coordinate of the AstroidalBall center.</span></div><div class='line' id='LC78'><span class="cm">     * @param a the factor of x.</span></div><div class='line' id='LC79'><span class="cm">     * @param b the factor of y.</span></div><div class='line' id='LC80'><span class="cm">     * @param c the factor of c.</span></div><div class='line' id='LC81'><span class="cm">     */</span></div><div class='line' id='LC82'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AstroidalBall</span><span class="p">(</span> <span class="k">const</span> <span class="kt">double</span> <span class="n">x0</span><span class="p">,</span> <span class="k">const</span> <span class="kt">double</span> <span class="n">y0</span><span class="p">,</span> <span class="k">const</span> <span class="kt">double</span>  <span class="n">z0</span><span class="p">,</span> <span class="k">const</span> <span class="kt">double</span> <span class="n">a</span><span class="p">,</span> <span class="k">const</span> <span class="kt">double</span> <span class="n">b</span><span class="p">,</span> <span class="k">const</span> <span class="kt">double</span>  <span class="n">c</span><span class="p">);</span></div><div class='line' id='LC83'><br/></div><div class='line' id='LC84'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="cm">/**</span></div><div class='line' id='LC85'><span class="cm">     * Constructor. </span></div><div class='line' id='LC86'><span class="cm">     * @param aPoint the AstroidalBall center.</span></div><div class='line' id='LC87'><span class="cm">     * @param aFactors the AstroidalBall factors</span></div><div class='line' id='LC88'><span class="cm">     */</span></div><div class='line' id='LC89'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AstroidalBall</span><span class="p">(</span><span class="k">const</span> <span class="n">RealPoint</span> <span class="o">&amp;</span><span class="n">aPoint</span><span class="p">,</span> <span class="k">const</span> <span class="n">RealPoint</span> <span class="o">&amp;</span><span class="n">aFactors</span><span class="p">);</span></div><div class='line' id='LC90'><br/></div><div class='line' id='LC91'>&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC92'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="cm">/**</span></div><div class='line' id='LC93'><span class="cm">     * Constructor. </span></div><div class='line' id='LC94'><span class="cm">     * @param aPoint the AstroidalBall center.</span></div><div class='line' id='LC95'><span class="cm">     * @param aFactors the AstroidalBall factors</span></div><div class='line' id='LC96'><span class="cm">     */</span></div><div class='line' id='LC97'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="cm">/*</span></div><div class='line' id='LC98'><span class="cm">    AstroidalBall(const Point &amp;aPoint, const RealPoint &amp;aFactors);</span></div><div class='line' id='LC99'><span class="cm">*/</span></div><div class='line' id='LC100'>&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC101'>&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC102'>&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC103'>&nbsp;&nbsp;<span class="c1">// ------------- Implementation of &#39;StarShaped&#39; services ------------------</span></div><div class='line' id='LC104'>&nbsp;&nbsp;<span class="k">public</span><span class="o">:</span></div><div class='line' id='LC105'><br/></div><div class='line' id='LC106'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="cm">/**</span></div><div class='line' id='LC107'><span class="cm">     * @return the lower bound of the AstroidalBall.</span></div><div class='line' id='LC108'><span class="cm">     *</span></div><div class='line' id='LC109'><span class="cm">     */</span></div><div class='line' id='LC110'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">RealPoint</span> <span class="n">getLowerBound</span><span class="p">()</span> <span class="k">const</span></div><div class='line' id='LC111'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC112'><br/></div><div class='line' id='LC113'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">return</span> <span class="n">RealPoint</span><span class="p">(</span><span class="n">myCenter</span><span class="p">[</span><span class="mi">0</span><span class="p">]</span> <span class="o">-</span> <span class="n">myFactors</span><span class="p">[</span><span class="mi">0</span><span class="p">],</span> <span class="n">myCenter</span><span class="p">[</span><span class="mi">1</span><span class="p">]</span> <span class="o">-</span> <span class="n">myFactors</span><span class="p">[</span><span class="mi">1</span><span class="p">],</span> <span class="n">myCenter</span><span class="p">[</span><span class="mi">2</span><span class="p">]</span> <span class="o">-</span> <span class="n">myFactors</span><span class="p">[</span><span class="mi">2</span><span class="p">]</span> <span class="p">);</span></div><div class='line' id='LC114'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC115'><br/></div><div class='line' id='LC116'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="cm">/**</span></div><div class='line' id='LC117'><span class="cm">     * @return the upper bound of the AstroidalBall.</span></div><div class='line' id='LC118'><span class="cm">     *</span></div><div class='line' id='LC119'><span class="cm">     */</span></div><div class='line' id='LC120'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">RealPoint</span> <span class="n">getUpperBound</span><span class="p">()</span> <span class="k">const</span></div><div class='line' id='LC121'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC122'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">return</span> <span class="n">RealPoint</span><span class="p">(</span><span class="n">myCenter</span><span class="p">[</span><span class="mi">0</span><span class="p">]</span> <span class="o">+</span> <span class="n">myFactors</span><span class="p">[</span><span class="mi">0</span><span class="p">],</span> <span class="n">myCenter</span><span class="p">[</span><span class="mi">1</span><span class="p">]</span> <span class="o">+</span> <span class="n">myFactors</span><span class="p">[</span><span class="mi">1</span><span class="p">],</span> <span class="n">myCenter</span><span class="p">[</span><span class="mi">2</span><span class="p">]</span> <span class="o">+</span> <span class="n">myFactors</span><span class="p">[</span><span class="mi">2</span><span class="p">]);</span></div><div class='line' id='LC123'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC124'><br/></div><div class='line' id='LC125'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="cm">/**</span></div><div class='line' id='LC126'><span class="cm">     * @return the center of the AstroidalBall.</span></div><div class='line' id='LC127'><span class="cm">     */</span></div><div class='line' id='LC128'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">RealPoint</span> <span class="n">center</span><span class="p">()</span> <span class="k">const</span></div><div class='line' id='LC129'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC130'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">return</span> <span class="n">myCenter</span><span class="p">;</span></div><div class='line' id='LC131'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC132'>&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC133'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="cm">/**</span></div><div class='line' id='LC134'><span class="cm">     * @param p any point in the space.</span></div><div class='line' id='LC135'><span class="cm">     *</span></div><div class='line' id='LC136'><span class="cm">     * @return the angle parameters wich are respectivly between [-Pi/2,Pi/2) and [-Pi,Pi] corresponding to</span></div><div class='line' id='LC137'><span class="cm">     * this point for the shape.</span></div><div class='line' id='LC138'><span class="cm">     */</span></div><div class='line' id='LC139'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AngularCoordinates</span> <span class="n">parameter</span><span class="p">(</span> <span class="k">const</span> <span class="n">RealPoint</span> <span class="o">&amp;</span> <span class="n">p</span> <span class="p">)</span> <span class="k">const</span><span class="p">;</span></div><div class='line' id='LC140'><br/></div><div class='line' id='LC141'><br/></div><div class='line' id='LC142'>&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC143'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="cm">/**</span></div><div class='line' id='LC144'><span class="cm">     * @param t is a couple of Teta &amp;&amp; Phi wich are respectivly between [-Pi/2,Pi/2) and [-Pi,Pi].</span></div><div class='line' id='LC145'><span class="cm">     *</span></div><div class='line' id='LC146'><span class="cm">     * @return the vector (x(t),y(t),z(t)) which is the position on the</span></div><div class='line' id='LC147'><span class="cm">     * shape boundary.</span></div><div class='line' id='LC148'><span class="cm">     */</span></div><div class='line' id='LC149'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">RealPoint</span> <span class="n">x</span><span class="p">(</span> <span class="k">const</span> <span class="n">AngularCoordinates</span> <span class="n">t</span> <span class="p">)</span> <span class="k">const</span><span class="p">;</span></div><div class='line' id='LC150'><br/></div><div class='line' id='LC151'><br/></div><div class='line' id='LC152'><br/></div><div class='line' id='LC153'><br/></div><div class='line' id='LC154'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="cm">/**</span></div><div class='line' id='LC155'><span class="cm">     * @param t is a couple of Teta &amp;&amp; Phi wich are respectivly between [-Pi/2,Pi/2) and [-Pi,Pi].</span></div><div class='line' id='LC156'><span class="cm">     *</span></div><div class='line' id='LC157'><span class="cm">     * @return the vector (gradf(M)).</span></div><div class='line' id='LC158'><span class="cm">     */</span></div><div class='line' id='LC159'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">virtual</span> <span class="n">RealPoint</span> <span class="n">gradient</span><span class="p">(</span> <span class="k">const</span> <span class="n">AngularCoordinates</span> <span class="n">t</span><span class="p">)</span> <span class="k">const</span> <span class="p">;</span></div><div class='line' id='LC160'><br/></div><div class='line' id='LC161'>&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC162'>&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC163'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="cm">/**</span></div><div class='line' id='LC164'><span class="cm">     * @param t is a couple of Teta &amp;&amp; Phi wich are respectivly between [-Pi/2,Pi/2) and [-Pi,Pi].</span></div><div class='line' id='LC165'><span class="cm">     *</span></div><div class='line' id='LC166'><span class="cm">     * @return the vector (rt(M)) wich is the partial derivative with respect to Teta.</span></div><div class='line' id='LC167'><span class="cm">     */</span></div><div class='line' id='LC168'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">virtual</span> <span class="n">RealPoint</span> <span class="n">rt</span><span class="p">(</span> <span class="k">const</span> <span class="n">AngularCoordinates</span> <span class="n">t</span><span class="p">)</span> <span class="k">const</span> <span class="p">;</span></div><div class='line' id='LC169'><br/></div><div class='line' id='LC170'><br/></div><div class='line' id='LC171'><br/></div><div class='line' id='LC172'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="cm">/**</span></div><div class='line' id='LC173'><span class="cm">     * @param t is a couple of Teta &amp;&amp; Phi wich are respectivly between [-Pi/2,Pi/2) and [-Pi,Pi].</span></div><div class='line' id='LC174'><span class="cm">     *</span></div><div class='line' id='LC175'><span class="cm">     * @return the vector (rp(M)) wich is the partial derivative with respect to Phi.</span></div><div class='line' id='LC176'><span class="cm">     */</span></div><div class='line' id='LC177'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">virtual</span> <span class="n">RealPoint</span> <span class="n">rp</span><span class="p">(</span> <span class="k">const</span> <span class="n">AngularCoordinates</span> <span class="n">t</span><span class="p">)</span> <span class="k">const</span> <span class="p">;</span></div><div class='line' id='LC178'><br/></div><div class='line' id='LC179'><br/></div><div class='line' id='LC180'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="cm">/**</span></div><div class='line' id='LC181'><span class="cm">     * @param t is a couple of Teta &amp;&amp; Phi wich are respectivly between [-Pi/2,Pi/2) and [-Pi,Pi].</span></div><div class='line' id='LC182'><span class="cm">     *</span></div><div class='line' id='LC183'><span class="cm">     * @return the vector (rtt(M)) wich is second the partial derivative with respect to Teta(twice).</span></div><div class='line' id='LC184'><span class="cm">     */</span></div><div class='line' id='LC185'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">virtual</span> <span class="n">RealPoint</span> <span class="n">rtt</span><span class="p">(</span> <span class="k">const</span> <span class="n">AngularCoordinates</span> <span class="n">t</span><span class="p">)</span> <span class="k">const</span> <span class="p">;</span></div><div class='line' id='LC186'><br/></div><div class='line' id='LC187'><br/></div><div class='line' id='LC188'><br/></div><div class='line' id='LC189'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="cm">/**</span></div><div class='line' id='LC190'><span class="cm">     * @param t is a couple of Teta &amp;&amp; Phi wich are respectivly between [-Pi/2,Pi/2) and [-Pi,Pi].</span></div><div class='line' id='LC191'><span class="cm">     *</span></div><div class='line' id='LC192'><span class="cm">     * @return the vector (rpp(M)) wich is second the partial derivatif with respect to Phi(twice).</span></div><div class='line' id='LC193'><span class="cm">     */</span></div><div class='line' id='LC194'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">virtual</span> <span class="n">RealPoint</span> <span class="n">rpp</span><span class="p">(</span> <span class="k">const</span> <span class="n">AngularCoordinates</span> <span class="n">t</span><span class="p">)</span> <span class="k">const</span> <span class="p">;</span></div><div class='line' id='LC195'>&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC196'><br/></div><div class='line' id='LC197'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="cm">/**</span></div><div class='line' id='LC198'><span class="cm">     * @param t is a couple of Teta &amp;&amp; Phi wich are respectivly between [-Pi/2,Pi/2) and [-Pi,Pi].</span></div><div class='line' id='LC199'><span class="cm">     *</span></div><div class='line' id='LC200'><span class="cm">     * @return the vector (rtp(M)) wich is second the partial derivatif with respect to Teta then phi.</span></div><div class='line' id='LC201'><span class="cm">     */</span></div><div class='line' id='LC202'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">virtual</span> <span class="n">RealPoint</span> <span class="n">rtp</span><span class="p">(</span> <span class="k">const</span> <span class="n">AngularCoordinates</span> <span class="n">t</span><span class="p">)</span> <span class="k">const</span> <span class="p">;</span></div><div class='line' id='LC203'>&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC204'>&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC205'>&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC206'>&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC207'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// ------------------------- data ----------------------------</span></div><div class='line' id='LC208'>&nbsp;&nbsp;<span class="k">private</span><span class="o">:</span></div><div class='line' id='LC209'><br/></div><div class='line' id='LC210'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="cm">/**</span></div><div class='line' id='LC211'><span class="cm">     * Center of the AstroidalBall.</span></div><div class='line' id='LC212'><span class="cm">     */</span></div><div class='line' id='LC213'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">RealPoint</span> <span class="n">myCenter</span><span class="p">;</span></div><div class='line' id='LC214'>&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC215'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="cm">/**</span></div><div class='line' id='LC216'><span class="cm">     * factors of the AstroidalBall.</span></div><div class='line' id='LC217'><span class="cm">     */</span></div><div class='line' id='LC218'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">RealPoint</span> <span class="n">myFactors</span><span class="p">;</span></div><div class='line' id='LC219'><br/></div><div class='line' id='LC220'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// ----------------------- Interface --------------------------------------</span></div><div class='line' id='LC221'>&nbsp;&nbsp;<span class="k">public</span><span class="o">:</span></div><div class='line' id='LC222'><br/></div><div class='line' id='LC223'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="cm">/**</span></div><div class='line' id='LC224'><span class="cm">     * Writes/Displays the object on an output stream.</span></div><div class='line' id='LC225'><span class="cm">     * @param out the output stream where the object is written.</span></div><div class='line' id='LC226'><span class="cm">     */</span></div><div class='line' id='LC227'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">void</span> <span class="n">selfDisplay</span> <span class="p">(</span> <span class="n">std</span><span class="o">::</span><span class="n">ostream</span> <span class="o">&amp;</span> <span class="n">out</span> <span class="p">)</span> <span class="k">const</span><span class="p">;</span></div><div class='line' id='LC228'><br/></div><div class='line' id='LC229'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="cm">/**</span></div><div class='line' id='LC230'><span class="cm">     * Checks the validity/consistency of the object.</span></div><div class='line' id='LC231'><span class="cm">     * @return &#39;true&#39; if the object is valid, &#39;false&#39; otherwise.</span></div><div class='line' id='LC232'><span class="cm">     */</span></div><div class='line' id='LC233'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">bool</span> <span class="n">isValid</span><span class="p">()</span> <span class="k">const</span><span class="p">;</span></div><div class='line' id='LC234'><br/></div><div class='line' id='LC235'><br/></div><div class='line' id='LC236'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// ------------------------- Hidden services ------------------------------</span></div><div class='line' id='LC237'>&nbsp;&nbsp;<span class="k">protected</span><span class="o">:</span></div><div class='line' id='LC238'><br/></div><div class='line' id='LC239'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="cm">/**</span></div><div class='line' id='LC240'><span class="cm">     * Constructor.</span></div><div class='line' id='LC241'><span class="cm">     * Forbidden by default (protected to avoid g++ warnings).</span></div><div class='line' id='LC242'><span class="cm">     */</span></div><div class='line' id='LC243'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AstroidalBall</span><span class="p">();</span></div><div class='line' id='LC244'><br/></div><div class='line' id='LC245'>&nbsp;&nbsp;<span class="k">private</span><span class="o">:</span></div><div class='line' id='LC246'><br/></div><div class='line' id='LC247'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="cm">/**</span></div><div class='line' id='LC248'><span class="cm">     * Copy constructor.</span></div><div class='line' id='LC249'><span class="cm">     * @param other the object to clone.</span></div><div class='line' id='LC250'><span class="cm">     * Forbidden by default.</span></div><div class='line' id='LC251'><span class="cm">     */</span></div><div class='line' id='LC252'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">//  AstroidalBall ( const AstroidalBall &amp; other );</span></div><div class='line' id='LC253'><br/></div><div class='line' id='LC254'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="cm">/**</span></div><div class='line' id='LC255'><span class="cm">     * Assignment.</span></div><div class='line' id='LC256'><span class="cm">     * @param other the object to copy.</span></div><div class='line' id='LC257'><span class="cm">     * @return a reference on &#39;this&#39;.</span></div><div class='line' id='LC258'><span class="cm">     * Forbidden by default.</span></div><div class='line' id='LC259'><span class="cm">     */</span></div><div class='line' id='LC260'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AstroidalBall</span> <span class="o">&amp;</span> <span class="k">operator</span><span class="o">=</span> <span class="p">(</span> <span class="k">const</span> <span class="n">AstroidalBall</span> <span class="o">&amp;</span> <span class="n">other</span> <span class="p">);</span></div><div class='line' id='LC261'><br/></div><div class='line' id='LC262'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// ------------------------- Internals ------------------------------------</span></div><div class='line' id='LC263'>&nbsp;&nbsp;<span class="k">private</span><span class="o">:</span></div><div class='line' id='LC264'><br/></div><div class='line' id='LC265'>&nbsp;&nbsp;<span class="p">};</span> <span class="c1">// end of class AstroidalBall</span></div><div class='line' id='LC266'><br/></div><div class='line' id='LC267'><br/></div><div class='line' id='LC268'>&nbsp;&nbsp;<span class="cm">/**</span></div><div class='line' id='LC269'><span class="cm">   * Overloads &#39;operator&lt;&lt;&#39; for displaying objects of class &#39;AstroidalBall&#39;.</span></div><div class='line' id='LC270'><span class="cm">   * @param out the output stream where the object is written.</span></div><div class='line' id='LC271'><span class="cm">   * @param object the object of class &#39;AstroidalBall&#39; to write.</span></div><div class='line' id='LC272'><span class="cm">   * @return the output stream after the writing.</span></div><div class='line' id='LC273'><span class="cm">   */</span></div><div class='line' id='LC274'>&nbsp;&nbsp;<span class="k">template</span> <span class="o">&lt;</span><span class="k">typename</span> <span class="n">T</span><span class="o">&gt;</span></div><div class='line' id='LC275'>&nbsp;&nbsp;<span class="n">std</span><span class="o">::</span><span class="n">ostream</span><span class="o">&amp;</span></div><div class='line' id='LC276'>&nbsp;&nbsp;<span class="k">operator</span><span class="o">&lt;&lt;</span> <span class="p">(</span> <span class="n">std</span><span class="o">::</span><span class="n">ostream</span> <span class="o">&amp;</span> <span class="n">out</span><span class="p">,</span> <span class="k">const</span> <span class="n">AstroidalBall</span><span class="o">&lt;</span><span class="n">T</span><span class="o">&gt;</span> <span class="o">&amp;</span> <span class="n">object</span> <span class="p">);</span></div><div class='line' id='LC277'><br/></div><div class='line' id='LC278'><span class="p">}</span> <span class="c1">// namespace DGtal</span></div><div class='line' id='LC279'><br/></div><div class='line' id='LC280'><br/></div><div class='line' id='LC281'><span class="c1">///////////////////////////////////////////////////////////////////////////////</span></div><div class='line' id='LC282'><span class="c1">// Includes inline functions.</span></div><div class='line' id='LC283'><span class="cp">#include &quot;AstroidalBall.ih&quot;</span></div><div class='line' id='LC284'><br/></div><div class='line' id='LC285'><span class="c1">//                                                                           //</span></div><div class='line' id='LC286'><span class="c1">///////////////////////////////////////////////////////////////////////////////</span></div><div class='line' id='LC287'><br/></div><div class='line' id='LC288'><span class="cp">#endif </span><span class="c1">// !defined AstroidBall_h</span></div><div class='line' id='LC289'><br/></div><div class='line' id='LC290'><span class="cp">#undef AstroidalBall_RECURSES</span></div><div class='line' id='LC291'><span class="cp">#endif </span><span class="c1">// else defined(AstroidalBall_RECURSES)</span></div></pre></div>
          </td>
        </tr>
      </table>
  </div>

          </div>
        </div>
      </div>
    </div>

  </div>

<div class="frame frame-loading large-loading-area" style="display:none;" data-tree-list-url="/AnisB/DGtal/tree-list/89108d7230c1510e0faf111f9128bb11a97b89e7" data-blob-url-prefix="/AnisB/DGtal/blob/89108d7230c1510e0faf111f9128bb11a97b89e7">
  <img src="https://a248.e.akamai.net/assets.github.com/images/spinners/octocat-spinner-64.gif?1338942895" height="64" width="64">
</div>

      </div>
      <div class="context-overlay"></div>
    </div>

      <div id="footer-push"></div><!-- hack for sticky footer -->
    </div><!-- end of wrapper - hack for sticky footer -->

      <!-- footer -->
      <div id="footer" >
        
  <div class="upper_footer">
     <div class="container clearfix">

       <!--[if IE]><h4 id="blacktocat_ie">GitHub Links</h4><![endif]-->
       <![if !IE]><h4 id="blacktocat">GitHub Links</h4><![endif]>

       <ul class="footer_nav">
         <h4>GitHub</h4>
         <li><a href="https://github.com/about">About</a></li>
         <li><a href="https://github.com/blog">Blog</a></li>
         <li><a href="https://github.com/features">Features</a></li>
         <li><a href="https://github.com/contact">Contact &amp; Support</a></li>
         <li><a href="https://github.com/training">Training</a></li>
         <li><a href="http://enterprise.github.com/">GitHub Enterprise</a></li>
         <li><a href="http://status.github.com/">Site Status</a></li>
       </ul>

       <ul class="footer_nav">
         <h4>Tools</h4>
         <li><a href="http://get.gaug.es/">Gauges: Analyze web traffic</a></li>
         <li><a href="http://speakerdeck.com">Speaker Deck: Presentations</a></li>
         <li><a href="https://gist.github.com">Gist: Code snippets</a></li>
         <li><a href="http://mac.github.com/">GitHub for Mac</a></li>
         <li><a href="http://windows.github.com/">GitHub for Windows</a></li>
         <li><a href="http://mobile.github.com/">Issues for iPhone</a></li>
         <li><a href="http://jobs.github.com/">Job Board</a></li>
       </ul>

       <ul class="footer_nav">
         <h4>Extras</h4>
         <li><a href="http://shop.github.com/">GitHub Shop</a></li>
         <li><a href="http://octodex.github.com/">The Octodex</a></li>
       </ul>

       <ul class="footer_nav">
         <h4>Documentation</h4>
         <li><a href="http://help.github.com/">GitHub Help</a></li>
         <li><a href="http://developer.github.com/">Developer API</a></li>
         <li><a href="http://github.github.com/github-flavored-markdown/">GitHub Flavored Markdown</a></li>
         <li><a href="http://pages.github.com/">GitHub Pages</a></li>
       </ul>

     </div><!-- /.site -->
  </div><!-- /.upper_footer -->

<div class="lower_footer">
  <div class="container clearfix">
    <!--[if IE]><div id="legal_ie"><![endif]-->
    <![if !IE]><div id="legal"><![endif]>
      <ul>
          <li><a href="https://github.com/site/terms">Terms of Service</a></li>
          <li><a href="https://github.com/site/privacy">Privacy</a></li>
          <li><a href="https://github.com/security">Security</a></li>
      </ul>

      <p>&copy; 2012 <span title="0.29024s from fe16.rs.github.com">GitHub</span> Inc. All rights reserved.</p>
    </div><!-- /#legal or /#legal_ie-->

      <div class="sponsor">
        <a href="http://www.rackspace.com" class="logo">
          <img alt="Dedicated Server" height="36" src="https://a248.e.akamai.net/assets.github.com/images/modules/footer/rackspaces_logo.png?1338942895" width="38" />
        </a>
        Powered by the <a href="http://www.rackspace.com ">Dedicated
        Servers</a> and<br/> <a href="http://www.rackspacecloud.com">Cloud
        Computing</a> of Rackspace Hosting<span>&reg;</span>
      </div>
  </div><!-- /.site -->
</div><!-- /.lower_footer -->

      </div><!-- /#footer -->

    

<div id="keyboard_shortcuts_pane" class="instapaper_ignore readability-extra" style="display:none">
  <h2>Keyboard Shortcuts <small><a href="#" class="js-see-all-keyboard-shortcuts">(see all)</a></small></h2>

  <div class="columns threecols">
    <div class="column first">
      <h3>Site wide shortcuts</h3>
      <dl class="keyboard-mappings">
        <dt>s</dt>
        <dd>Focus site search</dd>
      </dl>
      <dl class="keyboard-mappings">
        <dt>?</dt>
        <dd>Bring up this help dialog</dd>
      </dl>
    </div><!-- /.column.first -->

    <div class="column middle" style='display:none'>
      <h3>Commit list</h3>
      <dl class="keyboard-mappings">
        <dt>j</dt>
        <dd>Move selection down</dd>
      </dl>
      <dl class="keyboard-mappings">
        <dt>k</dt>
        <dd>Move selection up</dd>
      </dl>
      <dl class="keyboard-mappings">
        <dt>c <em>or</em> o <em>or</em> enter</dt>
        <dd>Open commit</dd>
      </dl>
      <dl class="keyboard-mappings">
        <dt>y</dt>
        <dd>Expand URL to its canonical form</dd>
      </dl>
    </div><!-- /.column.first -->

    <div class="column last" style='display:none'>
      <h3>Pull request list</h3>
      <dl class="keyboard-mappings">
        <dt>j</dt>
        <dd>Move selection down</dd>
      </dl>
      <dl class="keyboard-mappings">
        <dt>k</dt>
        <dd>Move selection up</dd>
      </dl>
      <dl class="keyboard-mappings">
        <dt>o <em>or</em> enter</dt>
        <dd>Open issue</dd>
      </dl>
      <dl class="keyboard-mappings">
        <dt><span class="platform-mac">⌘</span><span class="platform-other">ctrl</span> <em>+</em> enter</dt>
        <dd>Submit comment</dd>
      </dl>
    </div><!-- /.columns.last -->

  </div><!-- /.columns.equacols -->

  <div style='display:none'>
    <div class="rule"></div>

    <h3>Issues</h3>

    <div class="columns threecols">
      <div class="column first">
        <dl class="keyboard-mappings">
          <dt>j</dt>
          <dd>Move selection down</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>k</dt>
          <dd>Move selection up</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>x</dt>
          <dd>Toggle selection</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>o <em>or</em> enter</dt>
          <dd>Open issue</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt><span class="platform-mac">⌘</span><span class="platform-other">ctrl</span> <em>+</em> enter</dt>
          <dd>Submit comment</dd>
        </dl>
      </div><!-- /.column.first -->
      <div class="column last">
        <dl class="keyboard-mappings">
          <dt>c</dt>
          <dd>Create issue</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>l</dt>
          <dd>Create label</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>i</dt>
          <dd>Back to inbox</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>u</dt>
          <dd>Back to issues</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>/</dt>
          <dd>Focus issues search</dd>
        </dl>
      </div>
    </div>
  </div>

  <div style='display:none'>
    <div class="rule"></div>

    <h3>Issues Dashboard</h3>

    <div class="columns threecols">
      <div class="column first">
        <dl class="keyboard-mappings">
          <dt>j</dt>
          <dd>Move selection down</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>k</dt>
          <dd>Move selection up</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>o <em>or</em> enter</dt>
          <dd>Open issue</dd>
        </dl>
      </div><!-- /.column.first -->
    </div>
  </div>

  <div style='display:none'>
    <div class="rule"></div>

    <h3>Network Graph</h3>
    <div class="columns equacols">
      <div class="column first">
        <dl class="keyboard-mappings">
          <dt><span class="badmono">←</span> <em>or</em> h</dt>
          <dd>Scroll left</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt><span class="badmono">→</span> <em>or</em> l</dt>
          <dd>Scroll right</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt><span class="badmono">↑</span> <em>or</em> k</dt>
          <dd>Scroll up</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt><span class="badmono">↓</span> <em>or</em> j</dt>
          <dd>Scroll down</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>t</dt>
          <dd>Toggle visibility of head labels</dd>
        </dl>
      </div><!-- /.column.first -->
      <div class="column last">
        <dl class="keyboard-mappings">
          <dt>shift <span class="badmono">←</span> <em>or</em> shift h</dt>
          <dd>Scroll all the way left</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>shift <span class="badmono">→</span> <em>or</em> shift l</dt>
          <dd>Scroll all the way right</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>shift <span class="badmono">↑</span> <em>or</em> shift k</dt>
          <dd>Scroll all the way up</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>shift <span class="badmono">↓</span> <em>or</em> shift j</dt>
          <dd>Scroll all the way down</dd>
        </dl>
      </div><!-- /.column.last -->
    </div>
  </div>

  <div >
    <div class="rule"></div>
    <div class="columns threecols">
      <div class="column first" >
        <h3>Source Code Browsing</h3>
        <dl class="keyboard-mappings">
          <dt>t</dt>
          <dd>Activates the file finder</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>l</dt>
          <dd>Jump to line</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>w</dt>
          <dd>Switch branch/tag</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>y</dt>
          <dd>Expand URL to its canonical form</dd>
        </dl>
      </div>
    </div>
  </div>

  <div style='display:none'>
    <div class="rule"></div>
    <div class="columns threecols">
      <div class="column first">
        <h3>Browsing Commits</h3>
        <dl class="keyboard-mappings">
          <dt><span class="platform-mac">⌘</span><span class="platform-other">ctrl</span> <em>+</em> enter</dt>
          <dd>Submit comment</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>escape</dt>
          <dd>Close form</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>p</dt>
          <dd>Parent commit</dd>
        </dl>
        <dl class="keyboard-mappings">
          <dt>o</dt>
          <dd>Other parent commit</dd>
        </dl>
      </div>
    </div>
  </div>
</div>

    <div id="markdown-help" class="instapaper_ignore readability-extra">
  <h2>Markdown Cheat Sheet</h2>

  <div class="cheatsheet-content">

  <div class="mod">
    <div class="col">
      <h3>Format Text</h3>
      <p>Headers</p>
      <pre>
# This is an &lt;h1&gt; tag
## This is an &lt;h2&gt; tag
###### This is an &lt;h6&gt; tag</pre>
     <p>Text styles</p>
     <pre>
*This text will be italic*
_This will also be italic_
**This text will be bold**
__This will also be bold__

*You **can** combine them*
</pre>
    </div>
    <div class="col">
      <h3>Lists</h3>
      <p>Unordered</p>
      <pre>
* Item 1
* Item 2
  * Item 2a
  * Item 2b</pre>
     <p>Ordered</p>
     <pre>
1. Item 1
2. Item 2
3. Item 3
   * Item 3a
   * Item 3b</pre>
    </div>
    <div class="col">
      <h3>Miscellaneous</h3>
      <p>Images</p>
      <pre>
![GitHub Logo](/images/logo.png)
Format: ![Alt Text](url)
</pre>
     <p>Links</p>
     <pre>
http://github.com - automatic!
[GitHub](http://github.com)</pre>
<p>Blockquotes</p>
     <pre>
As Kanye West said:

> We're living the future so
> the present is our past.
</pre>
    </div>
  </div>
  <div class="rule"></div>

  <h3>Code Examples in Markdown</h3>
  <div class="col">
      <p>Syntax highlighting with <a href="http://github.github.com/github-flavored-markdown/" title="GitHub Flavored Markdown" target="_blank">GFM</a></p>
      <pre>
```javascript
function fancyAlert(arg) {
  if(arg) {
    $.facebox({div:'#foo'})
  }
}
```</pre>
    </div>
    <div class="col">
      <p>Or, indent your code 4 spaces</p>
      <pre>
Here is a Python code example
without syntax highlighting:

    def foo:
      if not bar:
        return true</pre>
    </div>
    <div class="col">
      <p>Inline code for comments</p>
      <pre>
I think you should use an
`&lt;addr&gt;` element here instead.</pre>
    </div>
  </div>

  </div>
</div>


    <div id="ajax-error-message">
      <span class="mini-icon mini-icon-exclamation"></span>
      Something went wrong with that request. Please try again.
      <a href="#" class="ajax-error-dismiss">Dismiss</a>
    </div>

    <div id="logo-popup">
      <h2>Looking for the GitHub logo?</h2>
      <ul>
        <li>
          <h4>GitHub Logo</h4>
          <a href="http://github-media-downloads.s3.amazonaws.com/GitHub_Logos.zip"><img alt="Github_logo" src="https://a248.e.akamai.net/assets.github.com/images/modules/about_page/github_logo.png?1338942895" /></a>
          <a href="http://github-media-downloads.s3.amazonaws.com/GitHub_Logos.zip" class="minibutton btn-download download"><span class="icon"></span>Download</a>
        </li>
        <li>
          <h4>The Octocat</h4>
          <a href="http://github-media-downloads.s3.amazonaws.com/Octocats.zip"><img alt="Octocat" src="https://a248.e.akamai.net/assets.github.com/images/modules/about_page/octocat.png?1338942895" /></a>
          <a href="http://github-media-downloads.s3.amazonaws.com/Octocats.zip" class="minibutton btn-download download"><span class="icon"></span>Download</a>
        </li>
      </ul>
    </div>

    
    
    
    <span id='server_response_time' data-time='0.29299' data-host='fe16'></span>
  </body>
</html>

