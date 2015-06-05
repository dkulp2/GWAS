###
# Compass
###

# Change Compass configuration
# compass_config do |config|
#   config.output_style = :compact
# end

###
# Page options, layouts, aliases and proxies
###

# Per-page layout changes:
#
# With no layout
# page "/path/to/file.html", :layout => false
#
# With alternative layout
# page "/path/to/file.html", :layout => :otherlayout
#
# A path which all have the same layout
# with_layout :admin do
#   page "/admin/*"
# end

# Proxy pages (https://middlemanapp.com/advanced/dynamic_pages/)
# proxy "/this-page-has-no-template.html", "/template-file.html", :locals => {
#  :which_fake_page => "Rendering a fake page with a local variable" }

###
# Helpers
###

# Automatic image dimensions on image_tag helper
# activate :automatic_image_sizes

# Reload the browser automatically whenever files change
configure :development do
  activate :livereload
end

helpers do
  # hack to set the active menu item for each page.  inspired by https://forum.middlemanapp.com/t/active-navigation/584/9
    def nav_link(link_text, url, options = {})
      options[:class] ||= ""
      active = (url == current_page.url ? ' active' : '')
      options[:class] << active
      "<li class='#{active}'>#{link_to(link_text, url, options)}</li>"
    end
#   def some_helper
#     "Helping"
#   end
end

set :markdown_engine, :kramdown

set :css_dir, 'stylesheets'

set :js_dir, 'javascripts'

set :images_dir, 'images'

activate :syntax

# Build-specific configuration
configure :build do
  # For example, change the Compass output style for deployment
  # activate :minify_css

  # Minify Javascript on build
  # activate :minify_javascript

  # Enable cache buster
  # activate :asset_hash

  # Use relative URLs
  activate :relative_assets

  # Or use a different image path
  #set :http_prefix, "/~dkulp"
end

activate :deploy do |deploy|
  deploy.method = :rsync
  deploy.host   = 'stat-gen.org'
  deploy.path   = '/home/dizzorg/public_html/statgen'
  # Optional Settings
  # deploy.user  = 'tvaughan' # no default
  # deploy.port  = 5309 # ssh port, default: 22
  # deploy.clean = true # remove orphaned files on remote host, default: false
  # deploy.flags = '-rltgoDvzO --no-p --del' # add custom flags, default: -avz
end
